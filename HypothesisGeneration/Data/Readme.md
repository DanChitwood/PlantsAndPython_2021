# Readme for Data Classification files

These files typically include a column with all of the unique names from one column of the original database, the count of rows with each unique name, and one or more columns with our new categories/classification for each unique value. See below for descriptions of each file.

## [AllTissues](AllTissues.csv)

95% of the data were categorized into tissue types, vegetative/reproduction, and above/below. 5% are still NAs, some with missing data and some that didn't fit nicely into categories. More cleaning may be necessary depending on the hypothesis being tested.

Columns and metadata:

- RawName: All unique names of tissue types in the original dataset. Includes all capitalizations, spelling, etc.
- Counts: Numbeer of rows that correspond to the tissue type RawName
- Tissue: Simplification of RawNames into one of 23 tissue types
- VegetativeRepro: Classifies the tissues into vegetative, reproductive, root, hypotocyl (spelling error - should be hypocotyl), or whole plant
- AboveBelow: Classifies the tissues into above, below, wholeplant, or seed. For the Above/Below hypothesis, only tissues with the values "Above" or "Below" should be used
- Debateable: Very little information in this column. "Yes" if the classifier (WL) wasn't completely sure how to categorize the tissue type. NA otherwise.
