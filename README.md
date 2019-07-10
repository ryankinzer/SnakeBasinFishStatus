
# Snake Basin Fish Status

Snake Basin Steelhead and Chinook population status assessment indicators and metrics.

## Goal

The goal of this work is to provide data for the upcoming status assessment. This should include estimates of total escapement to various TRT populations in the Snake River Basin, as well as estimates of female escapement (for productivity estimates) and estimates of escapement by brood year.

## Methods

We are using the [STADEM](https://github.com/KevinSee/STADEM) package to estimate total wild escapement past Lower Granite dam. We are using [PITcleanr](https://github.com/KevinSee/PITcleanr) and [DABOM](https://github.com/KevinSee/DABOM) to estimate total escapement to various tributaries.

## Data

Life history data (e.g. age, sex, etc.) comes from the Lower Granite trap database, which is maintained by IDFG. The valid tags from that trap database are then sent through a [PTAGIS](https://www.ptagis.org/) query to return their complete detection history. 
