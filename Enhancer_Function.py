import numpy as np
import pandas as pd
import anndata as ad
def ATAC_scEnhancer(adata1,scEnhancer):
    '''
    Calculate the matrix of synthetic lethal gene

    Parameters
    ----------
    None

    Returns
    -------
    ATAC_scEnhancer_csv:pandas.DataFrame
        the matrix of scATAC-seq index screened by the scEnhancer Data
    '''
    atac_var = adata1.var
    atac_var.columns = ['chrom', 'chromStart', 'chromEnd']

    atac_var
    chr = {'chr1_', 'chr2_', 'chr3_', 'chr4_', 'chr5_', 'chr6_', 'chr7_', 'chr8_', 'chr9_', 'chr10_', 'chr11_',
           'chr12_', 'chr13_', 'chr14_', 'chr15_', 'chr16_', 'chr17_', 'chr18_', 'chr19_', 'chr20_', 'chr21_', 'chr22_',
           'chrX_', 'chrY_'}

    scEnhancer.columns = ['chr', 'chromStart', 'chromEnd', 'Data']

    for i in range(len(scEnhancer.chr)):
        if (scEnhancer.chr.iloc[i] == 'chr1'):
            scEnhancer.chr.iloc[i] = 'chr1_'
        if (scEnhancer.chr.iloc[i] == 'chr2'):
            scEnhancer.chr.iloc[i] = 'chr2_'
        if (scEnhancer.chr.iloc[i] == 'chr3'):
            scEnhancer.chr.iloc[i] = 'chr3_'
        if (scEnhancer.chr.iloc[i] == 'chr4'):
            scEnhancer.chr.iloc[i] = 'chr4_'
        if (scEnhancer.chr.iloc[i] == 'chr5'):
            scEnhancer.chr.iloc[i] = 'chr5_'
        if (scEnhancer.chr.iloc[i] == 'chr6'):
            scEnhancer.chr.iloc[i] = 'chr6_'
        if (scEnhancer.chr.iloc[i] == 'chr7'):
            scEnhancer.chr.iloc[i] = 'chr7_'
        if (scEnhancer.chr.iloc[i] == 'chr8'):
            scEnhancer.chr.iloc[i] = 'chr8_'
        if (scEnhancer.chr.iloc[i] == 'chr9'):
            scEnhancer.chr.iloc[i] = 'chr9_'
        if (scEnhancer.chr.iloc[i] == 'chr10'):
            scEnhancer.chr.iloc[i] = 'chr10_'
        if (scEnhancer.chr.iloc[i] == 'chr11'):
            scEnhancer.chr.iloc[i] = 'chr11_'
        if (scEnhancer.chr.iloc[i] == 'chr12'):
            scEnhancer.chr.iloc[i] = 'chr12_'
        if (scEnhancer.chr.iloc[i] == 'chr13'):
            scEnhancer.chr.iloc[i] = 'chr13_'
        if (scEnhancer.chr.iloc[i] == 'chr14'):
            scEnhancer.chr.iloc[i] = 'chr14_'
        if (scEnhancer.chr.iloc[i] == 'chr15'):
            scEnhancer.chr.iloc[i] = 'chr15_'
        if (scEnhancer.chr.iloc[i] == 'chr16'):
            scEnhancer.chr.iloc[i] = 'chr16_'
        if (scEnhancer.chr.iloc[i] == 'chr17'):
            scEnhancer.chr.iloc[i] = 'chr17_'
        if (scEnhancer.chr.iloc[i] == 'chr18'):
            scEnhancer.chr.iloc[i] = 'chr18_'
        if (scEnhancer.chr.iloc[i] == 'chr19'):
            scEnhancer.chr.iloc[i] = 'chr19_'
        if (scEnhancer.chr.iloc[i] == 'chr20'):
            scEnhancer.chr.iloc[i] = 'chr20_'
        if (scEnhancer.chr.iloc[i] == 'chr21'):
            scEnhancer.chr.iloc[i] = 'chr21_'
        if (scEnhancer.chr.iloc[i] == 'chr22'):
            scEnhancer.chr.iloc[i] = 'chr22_'
        if (scEnhancer.chr.iloc[i] == 'chrX'):
            scEnhancer.chr.iloc[i] = 'chrX_'
        if (scEnhancer.chr.iloc[i] == 'chrY'):
            scEnhancer.chr.iloc[i] = 'chrY_'

    atac_var.chromStart = atac_var.chromStart.astype('int')
    atac_var.chromEnd = atac_var.chromEnd.astype('int')

    ATAC = {}
    ENHANCER = {}
    for x in chr:
        # 染色体拆分
        test = pd.DataFrame(columns=['Length'])
        ATAC[x] = test
        ATAC[x] = atac_var[atac_var.index.str.contains(x)]
        ENHANCER[x] = test
        test = pd.DataFrame(columns=['Length'])
        ENHANCER[x] = scEnhancer.loc[scEnhancer['chr'] == x]

        # 比对函数
        ATAC[x]['overleaf_1'] = ''
        ATAC[x]['overleaf_2'] = ''
        ATAC[x]['scEnhancer'] = ''
        for i in range(0, len(ENHANCER[x].chromEnd)):
            left = 0
            right = len(ATAC[x].chromEnd) - 1
            target_Start = ENHANCER[x].chromStart.iloc[i]
            target_End = ENHANCER[x].chromEnd.iloc[i]

            n = 0
            while (n < 30):
                # 中位数表示 #
                midIndex = round(int((left + right) / 2))
                midValue_Start = ATAC[x].chromStart.iloc[midIndex]
                midValue_End = ATAC[x].chromEnd.iloc[midIndex]
                # 三次condition判断，判断feature是否与增强子数据库匹配
                if (target_Start > midValue_End):
                    left = midIndex
                    n = n + 1
                elif (target_End < midValue_Start):
                    right = midIndex
                    n = n + 1

                else:
                    Length = target_End - target_Start
                    ATAC[x].overleaf_1.iloc[midIndex] = (target_Start - midValue_End) / Length
                    ATAC[x].overleaf_2.iloc[midIndex] = (midValue_Start - target_End) / Length
                    ATAC[x].scEnhancer.iloc[midIndex] = str(ENHANCER[x].chr.iloc[i]) + str(
                        ENHANCER[x].chromStart.iloc[i]) + '_' + str(ENHANCER[x].chromEnd.iloc[i])
                    break

    # 拼接dataframe
    ATAC_Enhancer = pd.concat(
        [ATAC['chr1_'], ATAC['chr2_'], ATAC['chr3_'], ATAC['chr4_'], ATAC['chr5_'], ATAC['chr6_'], ATAC['chr7_'],
         ATAC['chr8_'], ATAC['chr9_'], ATAC['chr10_'], ATAC['chr11_'], ATAC['chr12_'], ATAC['chr13_'], ATAC['chr14_'],
         ATAC['chr15_'], ATAC['chr16_'], ATAC['chr17_'], ATAC['chr18_'], ATAC['chr19_'], ATAC['chr20_'], ATAC['chr21_'],
         ATAC['chr22_'], ATAC['chrX_'], ATAC['chrY_'], ], axis=0)

    # 以下函数是删除overleaf未赋值的feature
    for i in range(0, len(ATAC_Enhancer.chromEnd)):
        if (ATAC_Enhancer.overleaf_1.iloc[i] == ''):
            ATAC_Enhancer.overleaf_1.iloc[i] = 'NaN'
            ATAC_Enhancer.overleaf_2.iloc[i] = 'NaN'
    ATAC_Enhancer

    # 保存索引，重置索引
    ATAC_Enhancer['index'] = ATAC_Enhancer.index
    ATAC_Enhancer.reset_index()

    # 比对，筛选
    ATAC_scEnhancer = pd.DataFrame(columns=['index', 'chrom', 'chromStart', 'chromEnd', 'scEnhancer'])
    n = 0
    for i in range(0, len(ATAC_Enhancer.chromEnd)):
        if (ATAC_Enhancer.overleaf_1.iloc[i] == 'NaN'):
            ATAC_Enhancer.overleaf_1.iloc[i] = 'NaN'
        else:
            ATAC_scEnhancer.loc[n] = ATAC_Enhancer.iloc[i]
            n = n + 1

    ATAC_scEnhancer
    ATAC_scEnhancer.index = ATAC_scEnhancer['index'].tolist()
    ATAC_scEnhancer.drop('index', axis=1, inplace=True)
    ATAC_scEnhancer_csv = ATAC_scEnhancer
    ATAC_scEnhancer_csv.to_csv('ATAC_scEnhancer.csv')
    return ATAC_scEnhancer_csv