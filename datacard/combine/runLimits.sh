DATACARDS=$1
SUBDIR=$(basename $DATACARDS)
if [[ "$DATACARDS" == "" ]]; then
    echo "ERROR: no datacard directory provided"
    exit 1
fi
mkdir -p workspaces/$SUBDIR
mkdir -p results/$SUBDIR
total=$(ls $DATACARDS/*.dat | wc -l)
counter=1
for datacard in $DATACARDS/*.dat; do
    name=$(basename ${datacard%%.dat})
    workspace=workspaces/$SUBDIR/${name}_workspace.root
    result=results/$SUBDIR/${name}_result.root
    log=results/$SUBDIR/${name}_logs.txt
    # Make workspace
    text2workspace.py $datacard -o $workspace &
done
echo "Waiting for workspaces to be created..."
wait


for datacard in $DATACARDS/*.dat; do
    name=$(basename ${datacard%%.dat})
    workspace=workspaces/$SUBDIR/${name}_workspace.root
    result=results/$SUBDIR/${name}_result.root
    log=results/$SUBDIR/${name}_logs.txt
    # Run limit
    combine -m 125 -M AsymptoticLimits -d ${workspace} -n .${name} --run blind &> $log &
    # combine -M Significance -m 125 -t -1 --expectSignal=1 -d ${workspace} -n .${name} &> $log &
done

echo "Waiting for limits to be calculated..."
wait

for datacard in $DATACARDS/*.dat; do
    name=$(basename ${datacard%%.dat})
    workspace=workspaces/$SUBDIR/${name}_workspace.root
    result=results/$SUBDIR/${name}_result.root
    log=results/$SUBDIR/${name}_logs.txt
    # Move result
    mv higgsCombine.${name}.AsymptoticLimits.mH125.root $result
    # mv higgsCombine.${name}.Significance.mH125.root $result
done

echo "Done!"