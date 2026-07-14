% columns:
%   1 - directory
%   2 - sex
%   3 - drug
%   4 - genetics (WT, AQ...)
%   5 - switch between matlab and Yue's model, using Yue's model as default
%   (now obsolete b/c we always use Yue's model)
%   6 - use FIE vs DLC diam (1 for DLC)

% from 8/21/24 email: C = CGRP
% from 11/3/25 email: P = PACAP

clear data_list
data_list{1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/OCT.2023/F15ANO-14102023'; data_list{end,2}='F'; data_list{end,3}='NO';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/OCT.2023/F16BNO-18102023'; data_list{end,2}='F'; data_list{end,3}='NO';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F28CNO3h-01122023'; data_list{end,2}='F'; data_list{end,3}='NO';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/OCT.2023/M21NO-24102023'; data_list{end,2}='M'; data_list{end,3}='NO';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/OCT.2023/M22NO-24102023'; data_list{end,2}='M'; data_list{end,3}='NO';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/OCT.2023/M24NO-28102023'; data_list{end,2}='M'; data_list{end,3}='NO';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/OCT.2023/F15BC-19102023'; data_list{end,2}='F'; data_list{end,3}='C';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F30C-12112023'; data_list{end,2}='F'; data_list{end,3}='C';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F26C-07112023'; data_list{end,2}='F'; data_list{end,3}='C';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M23BC-06112023'; data_list{end,2}='M'; data_list{end,3}='C';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M24BC-06112023'; data_list{end,2}='M'; data_list{end,3}='C';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M20BC-06112023'; data_list{end,2}='M'; data_list{end,3}='C';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M22BC-06112023'; data_list{end,2}='M'; data_list{end,3}='C';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F28AL-12112023'; data_list{end,2}='F'; data_list{end,3}='L';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/OCT.2023/M22AL-31102023'; data_list{end,2}='M'; data_list{end,3}='L';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/OCT.2023/M21AL-31102023'; data_list{end,2}='M'; data_list{end,3}='L';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M24AL-02112023'; data_list{end,2}='M'; data_list{end,3}='L';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M23AL-02112023'; data_list{end,2}='M'; data_list{end,3}='L';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/OCT.2023/M20AL-30102023'; data_list{end,2}='M'; data_list{end,3}='L';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F28BP-20112023'; data_list{end,2}='F'; data_list{end,3}='P';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M22CP-13112023'; data_list{end,2}='M'; data_list{end,3}='P';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M23CP-13112023'; data_list{end,2}='M'; data_list{end,3}='P';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M23CP-13112023'; data_list{end,2}='M'; data_list{end,3}='P';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F32M-12112023'; data_list{end,2}='F'; data_list{end,3}='M';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F31M-12112023'; data_list{end,2}='F'; data_list{end,3}='M';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F28M-07112023'; data_list{end,2}='F'; data_list{end,3}='M';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M35M-14112023'; data_list{end,2}='M'; data_list{end,3}='M';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M36M-14112023'; data_list{end,2}='M'; data_list{end,3}='M';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F31AS-22112023'; data_list{end,2}='F'; data_list{end,3}='S';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F32AS-22112023'; data_list{end,2}='F'; data_list{end,3}='S';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M22DS-27112023'; data_list{end,2}='M'; data_list{end,3}='S';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F33AA-23112023'; data_list{end,2}='F'; data_list{end,3}='A';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F28BA-27112023'; data_list{end,2}='F'; data_list{end,3}='A';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/F32BA-27112023'; data_list{end,2}='F'; data_list{end,3}='A';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M36AA-23112023'; data_list{end,2}='M'; data_list{end,3}='A';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M38A-23112023'; data_list{end,2}='M'; data_list{end,3}='A';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/DEC-2023/F32CSS3h-20122023'; data_list{end,2}='F'; data_list{end,3}='SS';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/FEB/AQF62BASELINE-E-LIGHT-130224'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/FEB/AQF60BASELINE-09022024'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/FEB/AQF64BASELINE-E-LIGHT-130224'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='AQ';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/WTF32A-200124'; data_list{end,2}='F'; data_list{end,3}='A'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/WTF42A-210124'; data_list{end,2}='F'; data_list{end,3}='A'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/AQF60A-190124'; data_list{end,2}='F'; data_list{end,3}='A'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/AQF62A-190124'; data_list{end,2}='F'; data_list{end,3}='A'; data_list{end,4}='AQ';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/WTF32KX-230124'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/WTF44KX-230124'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/AQF64KX-220124'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/AQF60KX-220124'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/AQF62KX-220124'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='AQ';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/FEB/AQF64PPA-020224'; data_list{end,2}='F'; data_list{end,3}='PPA'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/FEB/AQF62PPA-020224'; data_list{end,2}='F'; data_list{end,3}='PPA'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/FEB/AQF60PPA-050224'; data_list{end,2}='F'; data_list{end,3}='PPA'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/FEB/AQF62PPA-020224'; data_list{end,2}='F'; data_list{end,3}='PPA'; data_list{end,4}='AQ';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/AQF62NO-260124'; data_list{end,2}='F'; data_list{end,3}='NO'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/AQF60NO-290124'; data_list{end,2}='F'; data_list{end,3}='NO'; data_list{end,4}='AQ';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/FEB/AQF62C-080224'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/FEB/AQF64C-080224'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/AQM93C-300124'; data_list{end,2}='M'; data_list{end,3}='C'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/FEB/AQF60L-200224'; data_list{end,2}='F'; data_list{end,3}='L'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/FEB/AQF62P-260224'; data_list{end,2}='F'; data_list{end,3}='P'; data_list{end,4}='AQ';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/MARCH/M55PINP-11H-120324'; data_list{end,2}='M'; data_list{end,3}='PINP'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/MARCH/M53PINP-5H-120324'; data_list{end,2}='M'; data_list{end,3}='PINP'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/MARCH/M54pinp-8h-110324'; data_list{end,2}='M'; data_list{end,3}='PINP'; data_list{end,4}='WT';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/MAY 2024 AQPO4/AQM39-BASELINE-070524'; data_list{end,2}='M'; data_list{end,3}='BASELINE'; data_list{end,4}='AQ';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/MAY 2024 AQPO4/AQF412-BASELINE-160524'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/MAY 2024 AQPO4/AQF408-BASELINE-160524'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/AQM439-NO-070624'; data_list{end,2}='M'; data_list{end,3}='NO'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/MAY 2024 AQPO4/AQM439-KX-21052024'; data_list{end,2}='M'; data_list{end,3}='KX'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/MAY 2024 AQPO4/AQM443-BASELINE-020524'; data_list{end,2}='M'; data_list{end,3}='BASELINE'; data_list{end,4}='AQ';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/AQM447-KX-030624'; data_list{end,2}='M'; data_list{end,3}='KX'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/AQF406-NO-100624'; data_list{end,2}='F'; data_list{end,3}='NO'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/C57-WILDTYPE/C57-NO/F74-NO-210624'; data_list{end,2}='F'; data_list{end,3}='NO'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/C57-WILDTYPE/F71-BASELINE-120624'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/C57-WILDTYPE/C57-WILDTYPE MALE'; data_list{end,2}='M'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/C57-WILDTYPE/F70-BASELINE-130624'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/C57-WILDTYPE/C57-CGRP/F74-CGRP-250624'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/C57-WILDTYPE/F72-BASELINE-120624'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/AQM471-KX-110624'; data_list{end,2}='M'; data_list{end,3}='KX'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/AQF408-KX-040624'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/AQF406-KX-040624'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/AQM467-KX-110624'; data_list{end,2}='M'; data_list{end,3}='KX'; data_list{end,4}='AQ';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/C57-WILDTYPE/C57-CGRP/F70-CGRP-240624'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='WT';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/F80-BASELINE-030724'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/F87-NO-050724'; data_list{end,2}='F'; data_list{end,3}='NO'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/M92-DCLN-230724'; data_list{end,2}='M'; data_list{end,3}='DCLN'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/M95-DCLN-230724'; data_list{end,2}='M'; data_list{end,3}='DCLN'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/DCLN/M91-DCLF-060724'; data_list{end,2}='M'; data_list{end,3}='DCLF'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/DCLN/M92-DCLF-060724'; data_list{end,2}='M'; data_list{end,3}='DCLF'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/DCLN/M95-DCLF-060724'; data_list{end,2}='M'; data_list{end,3}='DCLF'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/F86-NO-050724'; data_list{end,2}='F'; data_list{end,3}='NO'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/F82-BASELINE-030724'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/M78-BASELINE-010724'; data_list{end,2}='M'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/F84-NO-040724'; data_list{end,2}='F'; data_list{end,3}='NO'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/F81-NO-020724'; data_list{end,2}='F'; data_list{end,3}='NO'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/O/F83-O-baseline-230724'; data_list{end,2}='F'; data_list{end,3}='O-baseline'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/O/F84-O-baseline-230724'; data_list{end,2}='F'; data_list{end,3}='O-baseline'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/F81-CGRP-060724'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/F83-CGRP-050724'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/F79-CGRP-060724'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/M78-NO-020724'; data_list{end,2}='M'; data_list{end,3}='NO'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/F82-CGRP-060724'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/F79-BASELINE-030724'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/CGRP/O/F87-O-baseline-230724'; data_list{end,2}='F'; data_list{end,3}='O-baseline'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/MAY 2024 AQPO4/AQM443-KX-170524'; data_list{end,2}='M'; data_list{end,3}='KX'; data_list{end,4}='AQ';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/C57-WILDTYPE/F74-BASELINE-130624'; data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/July/M77-NO-020724'; data_list{end,2}='M'; data_list{end,3}='NO'; data_list{end,4}='WT';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/AUG/F102-CGRP-210824'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='WT';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/AUG/M98-BASELINE-150824'; data_list{end,2}='M'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/F109-DEX-16092024'; data_list{end,2}='F'; data_list{end,3}='DEX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/M104-DEX-18092024'; data_list{end,2}='M'; data_list{end,3}='DEX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/M105-DEX-17092024'; data_list{end,2}='M'; data_list{end,3}='DEX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/M106-DEX-17092024'; data_list{end,2}='M'; data_list{end,3}='DEX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/M108-DEX-17092024'; data_list{end,2}='M'; data_list{end,3}='DEX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Oct/F109-DEX200-01102024'; data_list{end,2}='F'; data_list{end,3}='DEX200'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Oct/M117-DEX200-02102024'; data_list{end,2}='M'; data_list{end,3}='DEX200'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Oct/M104-DEX-INFUSION-04102024'; data_list{end,2}='M'; data_list{end,3}='DEX-INFUSION'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Oct/M108-DEX-INFUSION-03102024'; data_list{end,2}='M'; data_list{end,3}='DEX-INFUSION'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Oct/M106-DEX-INFUSION-03102024'; data_list{end,2}='M'; data_list{end,3}='DEX-INFUSION'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Oct/M116-DEX200-01102024'; data_list{end,2}='M'; data_list{end,3}='DEX200'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Oct/F109-DEX-INFUSION-07102024'; data_list{end,2}='F'; data_list{end,3}='DEX-INFUSION'; data_list{end,4}='WT';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/M104-BASELINE-060924'; data_list{end,2}='M'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/M108-baseline-05092024'; data_list{end,2}='M'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/M108-NO-09102024'; data_list{end,2}='M'; data_list{end,3}='NO'; data_list{end,4}='WT';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/F87-O-CGRP-220824 (repeat)'; data_list{end,2}='F'; data_list{end,3}='O-CGRP'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/F109-CGRP-12092024'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/F110-PACAP-20092024'; data_list{end,2}='F'; data_list{end,3}='P'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Oct/F84-O-CGRP-220824'; data_list{end,2}='F'; data_list{end,3}='O-CGRP'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Oct/F84-O-NO-240724'; data_list{end,2}='F'; data_list{end,3}='O-NO'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Oct/F87-O-NO-240724'; data_list{end,2}='F'; data_list{end,3}='O-NO'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/M106-NO-11092024'; data_list{end,2}='M'; data_list{end,3}='NO'; data_list{end,4}='WT';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F121-Dex-15112024'; data_list{end,2}='F'; data_list{end,3}='Dex'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F133-DEX-19112024'; data_list{end,2}='F'; data_list{end,3}='Dex'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/M128-DEX-06112024'; data_list{end,2}='M'; data_list{end,3}='Dex'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F118-DEX200-08112024'; data_list{end,2}='F'; data_list{end,3}='DEX200'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F119-DEX200-07112024'; data_list{end,2}='F'; data_list{end,3}='DEX200'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F121-DEX200-07112024'; data_list{end,2}='F'; data_list{end,3}='DEX200'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F133-DEX200-14112024'; data_list{end,2}='F'; data_list{end,3}='DEX200'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F134-DEX200-13112024'; data_list{end,2}='F'; data_list{end,3}='DEX200'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F135-DEX200-14112024'; data_list{end,2}='F'; data_list{end,3}='DEX200'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/M122-DEX200-04112024'; data_list{end,2}='M'; data_list{end,3}='DEX200'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/M128-DEX200-25112024'; data_list{end,2}='M'; data_list{end,3}='DEX200'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F131-Dex-15112024'; data_list{end,2}='F'; data_list{end,3}='Dex'; data_list{end,4}='WT';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F131-KX-28112024'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F132-KX-29112024'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F134-KX-29112024'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F135-KX-29112024'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/M127-KX-28112024'; data_list{end,2}='M'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/M129-KX-28112024'; data_list{end,2}='M'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/F142-KX-13122024'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/F144-KX-13122024'; data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/M138-KX-12122024'; data_list{end,2}='M'; data_list{end,3}='KX'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/M139-KX-12122024'; data_list{end,2}='M'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/M140-KX-12122024'; data_list{end,2}='M'; data_list{end,3}='KX'; data_list{end,4}='WT';
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/M148-kx-23122024'; data_list{end,2}='M'; data_list{end,3}='KX'; data_list{end,4}='WT';

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F153-Dex20-03022025'; data_list{end,2}='F'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F154-Dex20-07022025'; data_list{end,2}='F'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F156-Dex20-10022025'; data_list{end,2}='F'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F157-Dex20-05022025'; data_list{end,2}='F'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/F169-dex20-05032025'; data_list{end,2}='F'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/M161-dex20-12022025'; data_list{end,2}='M'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=1;
% M162-dex20 omitted
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/M165-dex20-05032025'; data_list{end,2}='M'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=1;

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F153-Dex40-06022025'; data_list{end,2}='F'; data_list{end,3}='Dex40'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F154-dex40-19022025'; data_list{end,2}='F'; data_list{end,3}='Dex40'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F155-dex40-16022025'; data_list{end,2}='F'; data_list{end,3}='Dex40'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F156-Dex40-04022025'; data_list{end,2}='F'; data_list{end,3}='Dex40'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F157-dex40-16022025'; data_list{end,2}='F'; data_list{end,3}='Dex40'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/F168-dex40-10032025'; data_list{end,2}='F'; data_list{end,3}='Dex40'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/M163-dex40-24022025'; data_list{end,2}='M'; data_list{end,3}='Dex40'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/M163-Dex40-26022025'; data_list{end,2}='M'; data_list{end,3}='Dex40'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/M166-dex40-04032025'; data_list{end,2}='M'; data_list{end,3}='Dex40'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/M167-dex40-04032025'; data_list{end,2}='M'; data_list{end,3}='Dex40'; data_list{end,4}='WT'; data_list{end,5}=1; data_list{end,6} = 1;

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F153-Dex60-10022025'; data_list{end,2}='F'; data_list{end,3}='Dex60'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F154-Dex60-04022025'; data_list{end,2}='F'; data_list{end,3}='Dex60'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F155-dex60-20022025'; data_list{end,2}='F'; data_list{end,3}='Dex60'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F157-dex60-20022025'; data_list{end,2}='F'; data_list{end,3}='Dex60'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/M161-dex60-19022025'; data_list{end,2}='M'; data_list{end,3}='Dex60'; data_list{end,4}='WT'; data_list{end,5}=1; data_list{end,6} = 1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/M163-dex60-03032025'; data_list{end,2}='M'; data_list{end,3}='Dex60'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/M165-dex60-03032025'; data_list{end,2}='M'; data_list{end,3}='Dex60'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/M166-dex60-25022025'; data_list{end,2}='M'; data_list{end,3}='Dex60'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/M167-dex60-27022025'; data_list{end,2}='M'; data_list{end,3}='Dex60'; data_list{end,4}='WT'; data_list{end,5}=1;

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/F142-amylin-20122024'; data_list{end,2}='F'; data_list{end,3}='amylin'; data_list{end,4}='WT'; data_list{end,5}=1;

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/F145-VIP-27122024'; data_list{end,2}='F'; data_list{end,3}='VIP'; data_list{end,4}='WT'; data_list{end,5}=1;

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/F169-dex20-13032025'; data_list{end,2}='F'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/F169-dex40-10032025'; data_list{end,2}='F'; data_list{end,3}='Dex40'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/M170-dex20-12032025'; data_list{end,2}='M'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/M171-dex20-12032025'; data_list{end,2}='M'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/M171-dex60-14032025'; data_list{end,2}='M'; data_list{end,3}='Dex60'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/M170-dex60-17032025'; data_list{end,2}='M'; data_list{end,3}='Dex60'; data_list{end,4}='WT'; data_list{end,5}=1;

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/M173-dex200-18032025'; data_list{end,2}='M'; data_list{end,3}='DEX200'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/M172-dex200-18032025'; data_list{end,2}='M'; data_list{end,3}='DEX200'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/F145-amylin-20122024'; data_list{end,2}='F'; data_list{end,3}='amylin'; data_list{end,4}='WT'; data_list{end,5}=1;

data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/Saline 2025/M170-saline-31032025'; data_list{end,2}='M'; data_list{end,3}='saline'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/Saline 2025/M173-saline26032025'; data_list{end,2}='M'; data_list{end,3}='saline'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/Saline 2025/M174-saline-26032025'; data_list{end,2}='M'; data_list{end,3}='saline'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/April-2025/F177-saline-01042025'; data_list{end,2}='F'; data_list{end,3}='saline'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/April-2025/F178-saline-04042025'; data_list{end,2}='F'; data_list{end,3}='saline'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/April-2025/F175-saline-02042025'; data_list{end,2}='F'; data_list{end,3}='saline'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/April-2025/F176-saline-01042025'; data_list{end,2}='F'; data_list{end,3}='saline'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/F144-VIP-26122024'; data_list{end,2}='F'; data_list{end,3}='VIP'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/F144-adm-18122024'; data_list{end,2}='F'; data_list{end,3}='adm'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/April-2025/F180-saline-11042025'; data_list{end,2}='F'; data_list{end,3}='saline'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/April-2025/F179-saline-03042025'; data_list{end,2}='F'; data_list{end,3}='saline'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/March-2025/Saline 2025/M171-saline-31032025'; data_list{end,2}='M'; data_list{end,3}='saline'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/Dec/F143-VIP-26122024'; data_list{end,2}='F'; data_list{end,3}='VIP'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2025/Feb-2025/F143-amylin-23122024'; data_list{end,2}='F'; data_list{end,3}='amylin'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/sept/F101-PACAP-12092024'; data_list{end,2}='F'; data_list{end,3}='P'; data_list{end,4}='WT'; data_list{end,5}=1;
%data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/sept/F102-PACAP-13092024'; data_list{end,2}='F'; data_list{end,3}='P'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/sept/F109-PACAP-23092024'; data_list{end,2}='F'; data_list{end,3}='P'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/sept/M106-PACAP-23092024'; data_list{end,2}='M'; data_list{end,3}='P'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/sept/M108-PACAP-24092024'; data_list{end,2}='M'; data_list{end,3}='P'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/sept/M105-CGRP-24092024'; data_list{end,2}='M'; data_list{end,3}='C'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/AUG/F103-CGRP-210824'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/AUG/M98-CGRP-200424'; data_list{end,2}='M'; data_list{end,3}='C'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/AUG/M99-CGRP-200424'; data_list{end,2}='M'; data_list{end,3}='C'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/NOV/M128-PACAP-12112024'; data_list{end,2}='M'; data_list{end,3}='P'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/NOV/M130-PACAP-12112024'; data_list{end,2}='M'; data_list{end,3}='P'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/Oct/M122-Maxi-25102024'; data_list{end,2}='M'; data_list{end,3}='M'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/sept/M104-NO-10092024'; data_list{end,2}='M'; data_list{end,3}='NO'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/JUNE/C57-WILDTYPE/C57-CGRP/F71-CGRP-250624'; data_list{end,2}='F'; data_list{end,3}='C'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/sept/M108-NO-11092024'; data_list{end,2}='M'; data_list{end,3}='NO'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/NOV/M127-PACAP-06112024'; data_list{end,2}='M'; data_list{end,3}='P'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/NOV/F131-adm-22112024'; data_list{end,2}='F'; data_list{end,3}='A'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/NOV/F131-maxi-18112024'; data_list{end,2}='F'; data_list{end,3}='M'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/NOV/F132-adm-21112024'; data_list{end,2}='F'; data_list{end,3}='A'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/NOV/F133-adm-26112024'; data_list{end,2}='F'; data_list{end,3}='A'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/NOV/F134-MAXI-20112024'; data_list{end,2}='F'; data_list{end,3}='M'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/NOV/M127-ADM-26112024'; data_list{end,2}='F'; data_list{end,3}='A'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/Dec/M138-ADM-17122024'; data_list{end,2}='M'; data_list{end,3}='A'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/NOV/F135-adm-21112024'; data_list{end,2}='F'; data_list{end,3}='A'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/NOV/M122-MAXI-18112024'; data_list{end,2}='M'; data_list{end,3}='M'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/sept/F102-DEX200-30092024'; data_list{end,2}='F'; data_list{end,3}='DEX200'; data_list{end,4}='WT'; data_list{end,5}=1;
data_list{end+1,1}='/mnt/nas/202309_Hashmat/2024/sept/F103-DEX200-30092024'; data_list{end,2}='F'; data_list{end,3}='DEX200'; data_list{end,4}='WT'; data_list{end,5}=1;

%find a way to chage the directory
old_path = '/gpfs/fs3/archive/dkell12_lab/202309_Hashmat';
new_path = '/mnt/nas/202309_Hashmat';
for i = 1:size(data_list,1)
    data_list{i,1} = strrep(data_list{i,1}, old_path, new_path);
end
%please create new data set after this line


% run thru all sets and set 4th column to WT if it is blank
for i = 1 : length(data_list)
    if isempty(data_list{i, 4})
        data_list{i, 4} = 'WT';
    end
end

% 5th column
% add extra column to flag whether to use yue's code
% look at ds = 156
for i = 1 : length(data_list)
    if isempty(data_list{i, 5})
        data_list{i, 5} = 0;
    end
end

% 6th column
% add extra column to flag whether to use diam (0) or diam_dlc (1)
for i = 1 : length(data_list)
    if isempty(data_list{i, 6})
        data_list{i, 6} = 0;
    end
% exclude the data sets you don't want to plot
% exclude data sets where the frequency was too low or part of the recording was bad
end

% initially we were making a list of data sets to plot (ds2plt) to exclude
% the data sets we didn't want to plot. We changed to not including these
% data sets in data_list, but for reverse compatability, we are still
% saving the variable.
ds2plt=1:size(data_list,1);

% exclude bad ECG data sets
exclude_ECG_sets = {'/mnt/nas/202309_Hashmat/2024/FEB/AQF64PPA-020224'; % bad ECG data sets
                    '/mnt/nas/202309_Hashmat/2024/FEB/AQF62PPA-020224';
                    '/mnt/nas/202309_Hashmat/2024/MAY 2024 AQPO4/AQF412-BASELINE-160524';
                    '/mnt/nas/202309_Hashmat/2024/MAY 2024 AQPO4/AQF408-BASELINE-160524';
                    '/mnt/nas/202309_Hashmat/2024/AUG/M98-BASELINE-150824';
                    '/mnt/nas/202309_Hashmat/2025/March-2025/F169-dex20-05032025';
                    '/mnt/nas/202309_Hashmat/2025/March-2025/M167-dex40-04032025';
                    };

exclude_ECG_mask = ~ismember(data_list(ds2plt,1),exclude_ECG_sets); % good ecg data sets
good_ECG_ds = ds2plt(exclude_ECG_mask);

%% make a list of excluded datasets
% these are data sets that were either commented out or are in the
% pulsatility projects status spreadsheet but were never completely analyzed b/c of image quality
% issues.
clear data_list_ex

%data sets that were once included but later commented out:
data_list_ex{1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/MAY 2024 AQPO4/AQF406-BASELINE-070524';% data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='AQ';
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/AQF412-KX-060624'; %data_list{end,2}='F'; data_list{end,3}='KX'; data_list{end,4}='AQ';
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/AUG/M99-BASELINE-19082024'; %data_list{end,2}='M'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/JUNE/C57-WILDTYPE/F73-BASELINE-200624'; %data_list{end,2}='F'; data_list{end,3}='BASELINE'; data_list{end,4}='WT';
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/sept/F112-DEX-19092024'; %data_list{end,2}='F'; data_list{end,3}='DEX'; data_list{end,4}='WT'; % don't use
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2024/NOV/F132-Dex-19112024'; %data_list{end,2}='F'; data_list{end,3}='Dex'; data_list{end,4}='WT';
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/F155-Dex20-05022025'; %data_list{end,2}='F'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=1; data_list{end,6} = 1;
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/April-2025/F179-saline-03042025'; %data_list{end,2}='F'; data_list{end,3}='saline'; data_list{end,4}='WT'; data_list{end,5}=0;
data_list_ex{end+1,1}='M158-dex20';
data_list_ex{end+1,1}='M160-dex20-11022025';

% omitted becuase almost no SVM
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/2025/Feb-2025/M160-dex20-18022025'; %data_list{end,2}='M'; data_list{end,3}='Dex20'; data_list{end,4}='WT'; data_list{end,5}=0; data_list{end,6}=1;

% frequency too low for cardiac, but could possibly be good for SVM
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/OCT.2023/F12NO+ECG'; data_list{end,2}='F'; data_list{end,3}='NO'; % frequency too low
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/OCT.2023/F16NO+ECG'; data_list{end,2}='F'; data_list{end,3}='NO'; % frequency too low
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M39NO3h-30112023'; data_list{end,2}='M'; data_list{end,3}='NO'; '202309_Hashmat/Nov.2023/M39NO3h-30112023'; % frequency too low (this was the weird one where the frequency changed during the recording)

% registration issues
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/Nov.2023/M38NO3h-29112023'; data_list{end,2}='M'; data_list{end,3}='NO';
data_list_ex{end+1,1}='/gpfs/fs3/archive/dkell12_lab/202309_Hashmat/DEC-2023/M36BNO3H-19122023'; data_list{end,2}='M'; data_list{end,3}='NO';

% data sets that are on the spreasheet but never finished due to image
% quality reasons
data_list_ex{end+1,1}='F15CL-23102023';
data_list_ex{end+1,1}='AQF64A-190124';
data_list_ex{end+1,1} = 'F27L-07112023';
data_list_ex{end+1,1} = 'F70-NO-210624';
data_list_ex{end+1,1} = 'AQF64NO-290124';
data_list_ex{end+1,1} = 'WTF44A-210124';
data_list_ex{end+1,1} ='AQ64L-E-LIGHT-220224';
data_list_ex{end+1,1} ='WTF32PPA-070224';
data_list_ex{end+1,1} ='WTF42PPA-050224';
data_list_ex{end+1,1} ='WTF44PPA-070224-New(080324)';
data_list_ex{end+1,1} ='WTF44PPA-070224-old';
data_list_ex{end+1,1}='AQM435-BASELINE-020524';
data_list_ex{end+1,1}='AQM445-PART1&2-BASELINE-020524';
data_list_ex{end+1,1}='M127-PACAP-11112024';
data_list_ex{end+1,1}='M127-PACAP-11112024';
DL = {'M129-DEX200-05112024'
'M129-PACAP-11112024'
'M130-DEX200-05112024'
'M130-MAXI-20112024'
'M129-ADM-25112024'
'F118-LEV-24102024'
'F119-LEV-11102024'
'F119-LEV-24102024'
'F120-LEV-10102024'
'F120-Lev-23102024'
'F121-Lev-23102024'
'M104-LEV-08102024'
'M108-LEV-08102024'
'M122-LEV-09102024'
'M122-LEV-22102024'
'M126-CGRP-21102024'
'F103-pacap-16092024'
'F111-PACAP-19092024'
'M106-CGRP-25092024'
'M116-CGRP-26092024'
'M117-CGRP-27092024'
'M162-dex20-18022025'
'F175-Gaunfacine-3mg-08042025'
'F175-Guanfacine-03042025'
'F176-Gaunfacine-02042025'
'F176-Guanfacine3mg-07042025'
'F177-Gaunfacine-3mg-08042025'
'F177-Guanfacine-02042025'
'F178-Gaunfacine-0,5mg-09042025'
'F179-Gaunfacine-0,5mg-10042025'
'F179-Gaunfacine-5mg-08042025'
'F180-Gaunfacine-5mg-09042025'
'M170-gaunfacine-0,5mg-09042025'
'M171-gaunfacine-0,5mg-09042025'
'M173-Gaunfacine 3mg-03042025'
'M173-Gaunfacine-5mg-08042025'
'M174-Gaunfacine-5mg-08042025'
'M174-Guanfacine-03042025'
'F34AS-22112023'
'M20CP-13112023'
'M24CM-13112023'
'M39A-23112023'
'M20NO-24102023'
'M23NO-28102023'
'M158-dex20-11022025'
};
for i = 1:length(DL)
    data_list_ex{end+1,1} = DL{i};
end

% The data sets that could work but going to take some extra efforts
DL_left = {
    'AQF408-NO-100624'
    'F71-NO-210624'
    'F143-adm-19122024'
    'F73_old'
    };
for i = 1:length(DL_left)
    data_list_ex{end+1,1} = DL_left{i};
end
DL_uncertain = {
    'AQF60-18032024'
    'AQF62-18032024'
    'AQF64-18032024'
    };
for i = 1:length(DL_uncertain)
    data_list_ex{end+1,1} = DL_uncertain{i};
end
Not_DL = {
    'Lecia Screen shots Jan-2024'
    'Lecia Screen shots Feb-2024'
    'Microscale slide Test 270224 Leica'
    '16-03-2024 AQF Test'
    '18032024'
    'Data251118_105250' %they are the same data set, F102 PACAP, noisy
    'F102-PACAP-13092024'%they are the same data set
    'Hashmat pilot May 2025'
    'Data_250420_220955'
    'Tiff files for EVA 25-04-2025'
    'F87-O-NO?'
    };
for i = 1:length(Not_DL)
    data_list_ex{end+1,1} = Not_DL{i};
end

%%
for i=1:length(data_list)
    if contains(data_list{i,3},'dex40','ignorecase',1)
        disp(data_list{i})
    end
end