"""Legacy-data additive benchmark, reusing BiMSGP fold metrics without rerunning it."""
from pathlib import Path
import argparse
import hashlib
import json
import re
import resource
import sys
import time
import numpy as np
import torch
from experiment import solve_grid, write_json, write_tsv, log
from preprocessing.panels import group_folds
from preprocessing.legacy import load_legacy, fit_transform, transform


def fit_context(x, y, train, test, models, rhos, config, directory, tag):
    output = {}
    for model in models:
        start = time.perf_counter()
        mean, feature_scale = fit_transform(x, train, model == 'Ridge-standardized')
        n = len(train)
        k = torch.zeros((n, n), dtype=torch.float64, device=config['device'])
        cross = torch.zeros((len(test), n), dtype=torch.float64, device=config['device'])
        for begin in range(0, x.shape[1], 4096):
            columns = np.arange(begin, min(begin + 4096, x.shape[1]))
            z = torch.as_tensor(transform(x, train, columns, mean, feature_scale), device=k.device)
            v = torch.as_tensor(transform(x, test, columns, mean, feature_scale), device=k.device)
            k.addmm_(z, z.T); cross.addmm_(v, z.T)
            del z, v
        kernel_scale = float(torch.trace(k) / n)
        k /= kernel_scale; cross /= kernel_scale
        prediction, residual, iterations = solve_grid(k, cross, y[train], rhos, config)
        output[model] = dict(prediction=prediction, residual=residual, iterations=iterations,
                             kernel_scale=kernel_scale, seconds=time.perf_counter() - start)
        np.savez(directory / f'{tag}_{model}_preprocessing.npz', train=train, test=test,
                 mean=mean, feature_scale=feature_scale, kernel_scale=kernel_scale)
        log(f'{directory.name} {tag} {model}: {len(train)} train; '
            f'{time.perf_counter()-start:.1f}s; max residual {max(residual):.2g}')
        del k, cross
    return output


def run_fold(f, inputs, config, output):
    x, y, ids, folds, markers, identity = inputs
    d = output / f'fold_{f+1}'; d.mkdir(exist_ok=True)
    train = np.flatnonzero(folds != f); test = np.flatnonzero(folds == f)
    # Fixed inner assignment without labels; the historical outer folds stay unchanged.
    inner = group_folds(np.array(ids)[train], config['inner_folds'], config['seed'], f'legacy-inner:{f}')
    write_tsv(d / 'inner_splits.tsv', [dict(sample_id=ids[i], inner_fold=int(k)) for i, k in zip(train, inner)])
    rhos = list(config['rho_grid']); started = time.perf_counter()
    if torch.cuda.is_available():torch.cuda.reset_peak_memory_stats()
    while True:
        scores = {m: np.zeros(len(rhos)) for m in config['models']}
        for k in range(config['inner_folds']):
            ti, vi = train[inner != k], train[inner == k]
            result = fit_context(x, y, ti, vi, config['models'], rhos, config, d, f'inner_{k}')
            for model, record in result.items():scores[model] += ((record['prediction'] - y[vi, None]) ** 2).sum(0)
        choice = {m: int(np.argmin(v)) for m, v in scores.items()}
        low = any(i == 0 for i in choice.values()) and rhos[0] > config['rho_bounds'][0]
        high = any(i == len(rhos)-1 for i in choice.values()) and rhos[-1] < config['rho_bounds'][1]
        if not low and not high:break
        if low:rhos.insert(0, max(rhos[0]/10, config['rho_bounds'][0]))
        if high:rhos.append(min(rhos[-1]*10, config['rho_bounds'][1]))
        log(f'fold {f+1}: extending shared search to {rhos[0]}..{rhos[-1]}')
    outer_rhos = sorted({rhos[j] for j in choice.values()})
    result = fit_context(x, y, train, test, config['models'], outer_rhos, config, d, 'outer')
    rows, search = [], []
    for model, record in result.items():
        rho = rhos[choice[model]]; j = outer_rhos.index(rho)
        rows += [dict(sample_id=ids[i], fold=f+1, model=model, observed=float(y[i]),
                      prediction=float(pred), rho=rho) for i, pred in zip(test, record['prediction'][:, j])]
        search.append(dict(model=model, rho=rho, search_rhos=rhos,
            inner_mse=(scores[model]/len(train)).tolist(),
            boundary_selected=choice[model] in [0,len(rhos)-1], relative_residual=record['residual'][j]))
    write_tsv(d / 'predictions.tsv', rows); write_json(d / 'search.json', search)
    write_json(d / 'complete.json', dict(status='COMPLETED', fold=f+1, train_samples=len(train),
        test_samples=len(test), seconds=time.perf_counter()-started,
        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
        peak_cuda_bytes=torch.cuda.max_memory_allocated() if torch.cuda.is_available() else 0))
    log(f'fold {f+1} completed in {(time.perf_counter()-started)/60:.2f}min')


def metrics(y, pred):
    return dict(mse=float(np.mean((y-pred)**2)), rmse=float(np.sqrt(np.mean((y-pred)**2))),
        pcc=float(np.corrcoef(y,pred)[0,1]), mae=float(np.mean(np.abs(y-pred))),
        r2=float(1-np.sum((y-pred)**2)/np.sum((y-y.mean())**2)))


def historical_metrics(path):
    text=path.read_text(); out={}
    for key in ['PCC','RMSE','MAE']:
        match=re.search(rf'{key} per fold: \[(.*?)\]',text)
        line=re.sub(r'np\.float64\(([^)]+)\)',r'\1',match.group(1))
        values=np.array(json.loads('['+line+']'),dtype=float)
        if len(values)!=10:raise ValueError('Historical fold metrics incomplete')
        out[key.lower()]=values
    out['mse']=out['rmse']**2
    return out


def report(inputs, config, root, output):
    import csv
    _, y, ids, folds, _, _=inputs
    rows=[]; costs=[]
    for f in range(10):
        d=output/f'fold_{f+1}'
        costs.append(json.loads((d/'complete.json').read_text()))
        with (d/'predictions.tsv').open() as s:rows+=list(csv.DictReader(s,delimiter='\t'))
    write_tsv(output/'predictions.tsv',rows)
    pooled=[]; fold_metrics=[]
    for model in config['models']:
        values={r['sample_id']:float(r['prediction']) for r in rows if r['model']==model}
        if set(values)!=set(ids):raise ValueError('Incomplete predictions')
        pred=np.array([values[i] for i in ids])
        pooled.append(dict(model=model,**metrics(y,pred)))
        for f in range(10):
            mask=folds==f;fold_metrics.append(dict(model=model,fold=f+1,samples=int(mask.sum()),**metrics(y[mask],pred[mask])))
    write_tsv(output/'pooled_metrics.tsv',pooled)
    write_tsv(output/'fold_metrics.tsv',fold_metrics)
    historical={name:historical_metrics((root/config[key]).resolve()) for name,key in
        [('BiMSGP-published','bimsgp_results'),('GBLUP-published','upstream_gblup_results')]}
    weights=np.bincount(folds,minlength=10); comparisons=[]
    for model in config['models']:
        current=[r for r in fold_metrics if r['model']==model]
        baseline_mse=np.array([r['mse'] for r in current])
        for reference, old in historical.items():
            historical_pooled_mse=float(np.dot(weights,old['mse'])/weights.sum())
            comparisons.append(dict(model=model,reference=reference,
                baseline_pooled_mse=float(np.dot(weights,baseline_mse)/weights.sum()),
                reference_pooled_mse_reconstructed=historical_pooled_mse,
                baseline_mean_fold_pcc=float(np.mean([r['pcc'] for r in current])),
                reference_mean_fold_pcc=float(old['pcc'].mean()),
                improved_mse_folds=int(np.sum(baseline_mse<old['mse'])),
                relative_mse_gain=float(1-np.dot(weights,baseline_mse)/np.dot(weights,old['mse'])),
                scope='descriptive matched-fold summary; historical epoch selection differs; no paired individual errors'))
    write_json(output/'historical_comparisons.json',comparisons)
    write_json(output/'costs.json',costs)
    write_json(output/'complete.json',dict(status='COMPLETED',samples=len(ids),markers=42195,
        outer_folds=10,baseline_encodings=config['models'],bimsgp_retrained=False,bimsgp_inference_run=False,
        official_sealed_label_file_read=False,scope=config['historical_cohort_scope']))
    lines=['# BiMSGP 原 SoyBase 队列的加性基准比较','',
        '使用既有 15,878 个材料、42,195 个原列序 SNP、protein 标签和原十折分配。BiMSGP 未重训、未运行推理，只复用公布的逐折指标。当前官方 VCF 密度实验保持独立。','',
        '| 本轮模型 | pooled MSE | RMSE | PCC | MAE |','| --- | ---: | ---: | ---: | ---: |']
    lines += [f'| {r["model"]} | {r["mse"]:.6f} | {r["rmse"]:.6f} | {r["pcc"]:.6f} | {r["mae"]:.6f} |' for r in pooled]
    lines += ['', 'GBLUP-centered 与 Ridge-standardized 属于同一加性岭回归模型类，区别是训练侧 SNP 方差缩放，不计作两个独立算法。两者均以内层 pooled MSE 选参，测试折不参与选参。', '',
        'BiMSGP 原成绩按测试折 PCC 选 epoch；本轮加性基准采用训练侧内层选参。数据和外折一致并不消除选模协议差异。历史汇总指标可做逐折描述性对照，不能恢复 pooled PCC、逐材料配对误差或重新构造无测试泄漏的 BiMSGP 结果。', '',
        '旧队列与官方 VCF 的封存名单有 3,139 个材料 ID 重合。本次按用户指定的完整历史队列使用旧标签，不读取官方封存标签文件；结果不能作为该官方队列的独立封存评价。', '',
        '完整预测、内外折统计、参数搜索及历史比较见本目录 TSV/JSON 文件。历史模型间精确性差异与成本比较须保留其协议和设备限制。']
    (output/'analysis.md').write_text('\n'.join(lines)+'\n')


def main():
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['profile','evaluate'])
    parser.add_argument('--config',default=str(Path(__file__).with_name('benchmark.json')))
    args=parser.parse_args();path=Path(args.config).resolve();root=path.parent;config=json.loads(path.read_text())
    # Reuse locally available parquet dependencies without installing model packages.
    sys.path.append(str(root.parent/'fprior'/'.deps'))
    from threadpoolctl import threadpool_limits
    threadpool_limits(config['blas_threads']);torch.set_num_threads(config['blas_threads'])
    inputs=load_legacy(config,root)
    output=(root/config['output']).resolve();output.mkdir(parents=True,exist_ok=True)
    sources=[Path(__file__),root/'preprocessing/legacy.py',root/'experiment.py',root/'preprocessing/panels.py']
    identity=dict(config=config,inputs=inputs[-1],source={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources})
    if (output/'identity.json').exists():
        if json.loads((output/'identity.json').read_text())!=identity:raise ValueError('Benchmark identity changed')
    else:write_json(output/'identity.json',identity)
    write_tsv(output/'marker_catalog.tsv',[dict(index=i,marker=m) for i,m in enumerate(inputs[4])])
    write_tsv(output/'outer_splits.tsv',[dict(sample_id=i,fold=int(f)+1) for i,f in zip(inputs[2],inputs[3])])
    log(f'Legacy data in RAM: {inputs[0].shape}, phenotype mean {inputs[1].mean():.4f}; BiMSGP not rerun')
    for f in ([0] if args.action=='profile' else range(10)):
        if (output/f'fold_{f+1}'/'complete.json').exists():continue
        run_fold(f,inputs,config,output)
    if args.action=='evaluate':report(inputs,config,root,output)


if __name__=='__main__':main()
