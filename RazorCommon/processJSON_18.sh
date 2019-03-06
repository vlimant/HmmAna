outDir=/mnt/hadoop/store/user/idutta/HmmAna/Data2018/JSON_processed_Dec042018
mkdir $outDir
for i in /mnt/hadoop/store/user/nlu/Hmm/ntuple/2018/Dec042018/Data/Prod_v1/Run2018*-Nano14Dec2018-v1_22Feb19_v2/*root;
do
    tempname="${i##*/}"
    filename="${tempname%.*}"
    outfile="${filename}_json.root"
    FWLiteGoodLumi Tools/python/loadJson_2018.py $i $outfile
    mv $outfile $outDir
done
