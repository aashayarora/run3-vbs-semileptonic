# Vector Boson Scattering - 1 lepton final state

## Preselection
Most corrections are applied using correctionlib, they can be found in `preselection/corrections/`. \
Preselection code is in `preselections/src/`, running `make` will compile and link all the files. The sample names are fed to the code using spec jsons in `preselections/input` which can be created using the `preselections/etc/make-config.py` script.

After creating the binary,
```
./bin/runAnalysis -i spec.json (--sig) (--bkg) (--data)
```

## MVA
`make` to compile and link the code in `mva/`
```
./bin/runMVA -i <input> -o <output>
```

## DNN
```
python3 python/infer.py configs/ABCDNet_simpleDisco_VBSVVH1lep_30.json --export --epoch 50
```

## Plotting
```
python3 plotting/make-plots.py
```