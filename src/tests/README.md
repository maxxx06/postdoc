# RIPTiDe

## How does riptide work ?

3 steps:

- First it reads the transcriptomic file. It takes in input a table like the one below then converts it in RPM (Read per million) :

| genes      | abundances |
| :---        |    :----:   |
|gene1     | nb1       | 
| gene2   | nb2        |
| ...   | ...        |


- Then, it assigns cefficients following the equation $\text{gw }=\frac{T_g}{T_max}$ where $\text{gw}$ represents reaction weights associated to gene g, ${T_g}$ represents the read-per-million normalised of the gene g and ${T_max}$ the  largest abundance of this distribution.

    -  check if GPR option is activated (Determines if GPR rules will be considered during coefficient assignment). 
        
        - **and** relationship: minimise transcrit abundance for the raction
        - **or** relationship: sum transcrit abundance for the reaction
        - **NONE**: Maximise transcrit abundance for the reaction

    - check if Additive option is activated (Pool transcription abundances for reactions with multiple contributing gene products)

        - **True**: sum transcrit abundance for the reaction
        - **False**: Maximise transcrit abundance for the reaction

- Finally, add constrains for **pruning** and **sampling** contrains.

    - *pruning*: multiply minimum coefficient of a reaction by its flux
    - *sampling*: multiply maximum coefficient of a reaction by its flux






