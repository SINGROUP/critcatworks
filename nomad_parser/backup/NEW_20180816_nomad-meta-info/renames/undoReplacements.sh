#!/bin/bash
# Removes all replacements to recover the old state
# call from this directory

reverse="--reverse"
#reverse=
if [ -e /labEnv3/bin/activate ]; then
    source /labEnv3/bin/activate
fi

cd ../meta_info_exploded/
echo -n "rename exploded... "
python ../renames/renamer.py $reverse */*.json
cd ../meta_info/meta_dictionary/
echo -n "rename compact ... "
python ../../renames/renamer.py $reverse *.meta_dictionary.json
cd ../nomad_meta_info
echo -n "rename v1      ... "
#python ../../renames/renamer.py $reverse *.nomadmetainfo.json
cd ../../..
sbt metaInfoTool/assembly
jar=meta-info-tool/target/scala-2.11/metaInfoTool-assembly-$(git describe --tags --always --dirty).jar
java -jar $jar reformat --v1 --base-dir=nomad-meta-info/meta_info/nomad_meta_info nomad-meta-info/meta_info/meta_dictionary/*.meta_dictionary.json
cd nomad-meta-info
find . -name \*.json | tar cf revertedMeta.tar -T -
git stash
git stash drop
git checkout v2_norename
tar xf revertedMeta.tar
git add -u .
rm revertedMeta.tar
git status
echo "moved to branch v2_norename, and ready to commit"
