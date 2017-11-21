#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -o grid_output
#$ -e grid_output
#$ -l h_vmem=5G
#$ -N ABSCAL_HERA_4POL

source activate HERA
# hand this script a SINGLE POL of files and it will get the others
ARGS=`pull_args.py $*`
PATH2PLIMPY='/home/saulkohn/githubs/plimpy'
PATH2HERACAL='/home/saulkohn/githubs/hera_cal'
path2abscal='/data4/paper/hera19_abscal/4polcal/2457548.46619.npz'
NOT_REAL_ANTS="0, 1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 13, 14, 15, 16, 17, 18, 19, 21, 23, 24, 25, 26, 27, 28, 29, 30, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, 52, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 66, 67, 68, 69, 70, 71, 73, 74, 75, 76, 77, 78, 79, 82, 83, 84, 85, 86, 87, 90, 91, 92, 93, 94, 95, 98, 99, 100, 101, 102, 103, 106, 107, 108, 109, 110, 111"

for f in ${ARGS}; do
    fQ=z`cut -d "z" -f 2 <<< "$f"`
    p=${fQ:18:2}
    fxx=${f/$p/xx}
    fxy=${f/$p/xy}
    fyx=${f/$p/yx}
    fyy=${f/$p/yy}
    # xx
    echo python ${PATH2PLIMPY}/reformat_and_apply_casa_gains.py -p xx --miriad=${fxx} --xy_npz=${path2abscal} --ex_ants="${NOT_REAL_ANTS}" --apply
    python ${PATH2PLIMPY}/reformat_and_apply_casa_gains.py -p xx --miriad=${fxx} --xy_npz=${path2abscal} --ex_ants="${NOT_REAL_ANTS}" --apply
    # xy
    echo python ${PATH2PLIMPY}/reformat_and_apply_casa_gains.py -p xy --miriad=${fxy} --xy_npz=${path2abscal} --ex_ants="${NOT_REAL_ANTS}" --apply
    python ${PATH2PLIMPY}/reformat_and_apply_casa_gains.py -p xy --miriad=${fxy} --xy_npz=${path2abscal} --ex_ants="${NOT_REAL_ANTS}" --apply
    # yx
    echo python ${PATH2PLIMPY}/reformat_and_apply_casa_gains.py -p yx --miriad=${fyx} --xy_npz=${path2abscal} --ex_ants="${NOT_REAL_ANTS}" --apply
    python ${PATH2PLIMPY}/reformat_and_apply_casa_gains.py -p yx --miriad=${fyx} --xy_npz=${path2abscal} --ex_ants="${NOT_REAL_ANTS}" --apply
    # yy
    echo python ${PATH2PLIMPY}/reformat_and_apply_casa_gains.py -p yy --miriad=${fyy} --xy_npz=${path2abscal} --ex_ants="${NOT_REAL_ANTS}" --apply
    python ${PATH2PLIMPY}/reformat_and_apply_casa_gains.py -p yy --miriad=${fyy} --xy_npz=${path2abscal} --ex_ants="${NOT_REAL_ANTS}" --apply
done

