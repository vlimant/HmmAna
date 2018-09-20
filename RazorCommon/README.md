1. In the Tools/bin/FWLiteGoodLumi.cc file, modify

   std::string RUN_BRANCH = "run";
   std::string LUMI_BRANCH = "lumi";

   To the appropriate branch names in your tree
2. In Tools/python/loadJSON.py file, add the correct name/location of your JSON input

3. In the processJSON.sh file, make sure to have the correct input and output directory names

