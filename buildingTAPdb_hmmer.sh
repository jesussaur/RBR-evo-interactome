# To concatenate hmm3 profiles 
# into 
# a local database file
cd Documents/Hmmer
for i in *.hmm3
do 
cat $i >> TAPdb;done


