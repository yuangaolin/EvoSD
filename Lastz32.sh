# #! /bin/bash
# number="$1"
# #b=$(number)
# #echo $b
# for ((i=1; i<=number; i++)); do
# 	/home/fanghuixuan/lastz-distrib/bin/lastz_32 AllSegsOfSCN_CRRemoved.fasta[multiple] AllSegsOfSCN_CRRemoved$i.fasta --format=general:cigar,score,name1,strand1,size1,zstart1,end1,name2,strand2,size2,zstart2,end2,identity,coverage,nmatch,nmismatch,ngap,cgap  > AllSegsOfSCN_CRRemoved$i.info &
# done
# wait
# echo end lastz alignments of putative SDs

#--notransition --step=20 --nogapped
#!/bin/bash

number="$1"

# 同时最多运行几个 LASTZ
MAX_JOBS=16

LASTZ=/home/fanghuixuan/lastz-distrib/bin/lastz_32

MASTER=AllSegsOfSCN_CRRemoved.fasta

FAILED=failed_lastz.txt

# 清空失败列表
> "$FAILED"


run_one() {

    i="$1"

    query="AllSegsOfSCN_CRRemoved${i}.fasta"

    out="AllSegsOfSCN_CRRemoved${i}.info"

    tmp="${out}.tmp"

    echo "START $i  $(date)"

    # 先写 tmp
    # 只有 LASTZ 正常结束后才替换正式 .info
    if "$LASTZ" \
        "${MASTER}[multiple]" \
        "$query" \
        --format=general:cigar,score,name1,strand1,size1,zstart1,end1,name2,strand2,size2,zstart2,end2,identity,coverage,nmatch,nmismatch,ngap,cgap \
        > "$tmp"
    then

        mv "$tmp" "$out"

        echo "OK    $i  $(date)"

    else

        status=$?

        rm -f "$tmp"

        echo "$i" >> "$FAILED"

        echo "FAILED $i exit=$status $(date)" >&2

    fi
}


for ((i=1; i<=number; i++))
do

    # 等待，直到当前后台任务 < MAX_JOBS
    while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]
    do
        sleep 2
    done

    run_one "$i" &

done

wait


echo
echo "========================================"

if [ -s "$FAILED" ]
then
    echo "Some LASTZ jobs FAILED:"
    cat "$FAILED"
else
    echo "All LASTZ jobs completed successfully."
fi

echo "========================================"

echo "end lastz alignments of putative SDs"
