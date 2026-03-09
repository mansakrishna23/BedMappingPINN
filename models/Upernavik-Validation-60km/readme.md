#### Validation Experiment: Coupling mass balance (MB) and stress balance (SB)
PINN (V3) trained, with different random seeds, for the Upernavik region using only a subset of ice thickness training data, spaced at ~60 km apart

Each subfolder contains the following files: 
- **history.json** : the PINN training history
- **history.png** : Figure of the PINN training history for different loss function componenets
- **model-700000.weights.h5** : PINN model weights saved after 700000 iterations
- **params.json** : Parameter file, can be used for setting up or re-loading the PINN
