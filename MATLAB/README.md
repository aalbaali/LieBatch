# What's in this repo?
Batch estimation on Lie groups. The data is assumed to be provided by the user and it has the following structure:
```
data_struct
    ├── sim
    |   ├── time
    |   └── freq    
    └── meas
        ├── prior
        |       ├── mean
        |       ├── cov
        |       └── tim
        ├── velocity
        |       ├── mean
        |       ├── cov
        |       ├── tim
        |       └── freq
        ├── gyro
        |       ├── mean
        |       ├── cov
        |       ├── tim
        |       └── freq
        └── gps
                ├── mean
                ├── cov
                ├── tim
                └── freq
```
### Notes about the data
- `cov` is a 3D array where the 3rd dimension is temporal (i.e., in time).
- `freq` is the data frequency.
- More data can be added to the struct but the code must then be modified.