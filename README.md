# OVERLAP-TOOLKIT

Utility tool(s) for protein-ligand modeling

## o-lap

Prune out overlapping atoms from mol2-formatted file
```
o-lap model.mol2
```

![o-lap modeling](https://media.springernature.com/lw685/springer-static/image/art%3A10.1186%2Fs13321-024-00857-6/MediaObjects/13321_2024_857_Figa_HTML.png)
### Usage

```
Usage: o-lap [options] model
Remove overlapping atoms from a model.

May use Markov Cluster Algorithm (MCL) tool for clustering.
Output is either a mol2 model, input for MCL, or atom types in MCL clusters.

Options:
  -h, --help             Displays help on commandline options.
  --help-all             Displays help including Qt specific options.
  -v, --version          Displays version information.
  --cutoffs <file/json>  JSON formatted cutoffs for atom types. Defaults are
                         are from file 'cutoffs.json'.
  -c, --cutoff <num>     Cutoff distance. Effective only when shorter than
                         default atomtype specific values (default: 1.1)
  --showcutoffs          Show JSON formatted cutoffs for atom types and exit.
  --similarjson <file>   JSON formatted atom types.
  --showsimilar          Show similar atom types and exit.
  -s, --similar          Cluster similar atom types. Types are similar, if in
                         same category. See '--showsimilar'
  --chargediff <num>     Charges must be within <num> to cluster (default:
                         0.2).
  --deletetypes <str>    Comma-separated list of atom types to discard
                         completely.
  --nib                  Convert atom types to (positive) N.3, (negative) O.3,
                         and (neutral) C.3/C.ar based on charge in NIB-like
                         manner.
  --nibthreshold <num>   Threshold of charge to bin atoms into N, C, O classes
                         (default: 0.2).
  --nibcharged           Keep all charges in model processed with nib option.
  --clustermin <int>     Minimum size of cluster to include (default: 1).
  --clusterminchr <int>  Minimum size of cluster for charged atoms. Atom is
                         charged, if abs(charge) exceeds nibthreshold. (default:
                         clustermin)
  --abcout               Create ABC-format input for MCL and exit.
  --mcl                  Create ABC-format input for MCL and run MCL.
  --mclI <num>           MCL main inflation value.
  --mclte <int>          MCL expansion thread number.
  --mapmcl <file>        Map mcl clusters to atoms. The <file> must be output
                         from MCL that corresponds to the model.
  --mcltype              Show types of clustered atoms.  Requires mapmcl.
  --prefix <str>         Prefix of the output molecule's name (default: model).

Arguments:
  model                  Mol2-file
```


Option `--cutoffs` requires a file or a value that is in JSON format. An example for latter:
```
o-lap --cutoffs '{"O.co2":2.2,"C.ar":1.1}' model.mol2
```
will use different cutoff for two atomtypes (and program's defaults for other types):
| Type | Cutoff |
| --- | ----------- |
| O.co2 | 2.2 |
| C.ar | 1.1 |

The default atom typing and cutoffs are in `INSTALL_PREFIX/share/SBL/o-lap/`.

## Dependencies

* [Qt 5](https://www.qt.io/): application and UI framework
* [MCL](https://micans.org/mcl/): Markov cluster algorithm (optional)

## License

GPLv3 (or later). See [LICENSE](https://github.com/jvlehtonen/overlap-toolkit/blob/main/LICENSE).

## How to build

Build requires C++ compiler that supports **C++17**, Qt 5, and CMake (version 3.10 or later)

### Build and install with CMake
```
cd overlap-toolkit
cmake -B build -DCMAKE_INSTALL_PREFIX=mypath
cmake --build build
cmake --install build
```

The `CMAKE_INSTALL_PREFIX` is `~/.local` by default.
Override that with `-DCMAKE_INSTALL_PREFIX=mypath` to direct installation into `mypath`.

If you did install to `mypath`, then `mypath/bin` must be on `PATH`.
For option `--mcl` the Markov cluster algorithm program `mcl` must be on `PATH`.


## How to cite Overlap Toolkit methods

1. For o-lap:
      Moyano-Gómez et al. 2024; J Cheminform 16, 97 (2024).
      doi: 10.1186/s13321-024-00857-6; https://doi.org/10.1186/s13321-024-00857-6
