# Parameter File for LBFlow

The parameter file is of [JSON] format.

[JSON]: http://www.json.org/ 

## Limits

 - Maximum length of parameter file: 8Kb.
 - Maximum number of characters in a string (such as file name): 255.
 - Gray scale range in data file: 0 ~ 255.

Find definition of these limits in filepar.h

## Keywords

Keywords of pure integer is considered as grayscale value.

### Grid definition

* **gridFile** - grid file name (_the path is "." by default_).
* **dimension** - grid dimension, in {x, y, z} order.
* **directions** - directions, i.e. 19 for D3Q19 and 9 for D2Q9.
* **latticeLen** - (m) side length of a unit lattice.

~~~~~~~~~~~~
"grid" : {
    "gridFile":     "Channel3D-24x24y100z.raw",
    "dimension":    ["24", "24", "100"],
    "directions":   "19",
    "latticeLen":   "1e-6"
},
~~~~~~~~~~~~

### Simulation control

* **startTime** - (ts) the start time step for the simulation.
                  (_For restart, startTime should be equal to the number of the .rst file, the previous output will not be overwritten._)
* **endTime** - (ts) the end time step for the simulation.
* **writeInterval** - (ts) the time step interval for writing output file.
* **reportInterval** - (ts) the time step interval for report.
* **restartInterval** - (ts) the time step interval for writing restart file.
* **maxError** - the error threshold for termination of the simulation.

~~~~~~~~~~~~
"control" : {
    "startTime":        "0",
    "endTime":          "1000",
    "writeInterval":    "100",
    "reportInterval":   "10",
    "restartInterval":  "500",
    "maxError":         "0.000001"
},
~~~~~~~~~~~~

### Output control

* **outPath** - the path and file name prefix for output.
* **reportDir** - direction for report.
* **reportSlice** - positions of report slice.

~~~~~~~~~~~~
"output" : {
    "outPath":      "./LBFlow-Result/Channel3D-24x24y100z",
    "reportDir":    "1",
    "reportSlice":  ["1", "98"]
},
~~~~~~~~~~~~

### Model

* **fluidModel** - case-insensitive name of the fluid model.
  - scsp: single component single phase model.
  - color: color gradient model for binary fluid.
  - sc: multiple component and multiple phase.
* **equilModel** - case-insensitive equilibrium function.
  - SRT: single rexlation time model.
  - MRT: multi-rexlation time model.
* **rexlationTime** - (ts) the relaxation time.
* **latticeSpeed** - (lu/ts) the lattice speed.
* **alpha0** - the free parameter determins the sound speed.
* **density** - (mu*lu^{-3}) the fluid density.
* **interfacialTension** - (mu/ts^2) the interfacial tension.
* **interfaceThickness** - this parameter controls the interface thickness.

~~~~~~~~~~~~
"model" : {
    "fluidModel":           "COLOR",
    "equilModel":           "MRT",
    "rexlationTime":        ["0.57", "1.215416667"],
    "latticeSpeed":         "1",
    "alpha0":               "0.33333333333333333333",
    "density":              ["1", "0.831683168"],
    "interfacialTension":   "0.085920139",
    "interfaceThickness":   "0.7"
},
~~~~~~~~~~~~

### Index

* **grayscale value** - the grayscale value corresponding with the luminance intensity in a 8-bit grid/image file (_grayscale must be number between 0 ~ 255_).
* **type**
  * "VD" - void;
  * "GV" - void with gravity settings;
  * "PB" - pressure boundary;
  * "VB" - velocity boundary;
  * "SD" - solid, apply contact angle of phase 1 by value.
           Periodic boundary condition is automatically applied.
           Use solid boundary to avoid periodicity.
* **value** - value of certain gird type (float).
* **dir** - direction for certain grid type.
  * Force direction follows DmQn scheme.
    * For D3Q19 1:x+, 2:x-, 3:y+, 4:y-, 5:z+, 6:z-;
    * for D2Q9 1:x+, 2:y+, 3:x-, 4: y-.
  * For velocity and pressure boundary, force direction must point to the inside, e.g. at east boundary, dir should point to the west.
  * For gravity, direction could be any direction.
* **phase** - phase index.

#### Required values

| Type | value | dir | phase |
|------|-------|-----|-------|
|  VD  |       |     |   Y   |
|  GV  |   Y   |  Y  |       |
|  PB  |   Y   |  Y  |   Y   |
|  VB  |   Y   |  Y  |   Y   |
|  SD  |   Y   |     |       |

~~~~~~~~~~~~
"index" : {
    "0" : {
        "type":  "VD",
        "phase": "0"
    },
    "50" : {
        "type":  "VB",
        "value": "0.00001",
        "dir":   "1",
        "phase": "0"
    },
    "75" : {
        "type":  "GV",
        "value": "0.001",
        "dir":   "1"
    },
    "200" : {
        "type":  "PB",
        "value": "1",
        "dir":   "2",
        "phase": "0"
    },
    "255" : {
        "type":  "SD",
        "value": "160"
    }
}
~~~~~~~~~~~~

