outDir=/eos/user/i/idutta/JSON_processed_20Sep18
mkdir $outDir
for i in /eos/cms/store/user/nlu/Hmm/ntuple/Data2017*root;
do
    tempname="${i##*/}"
    filename="${tempname%.*}"
    outfile="${filename}_json.root"
    FWLiteGoodLumi Tools/python/loadJson.py $i $outfile
    mv $outfile $outDir
done
