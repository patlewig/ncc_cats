# Input Formats and examples

## Input Format

The classification methods rely on the Distributed Structure Searchable Toxicity Substance Identifier DSSTOX_SID, the log10 of the octanol-water partition coefficient (Log Kow), Molecular Weight, Simplified Molecular Input Line Entry System (SMILES), Water Solubility,and a molecular object (RDkit.Mol) generated using the python library RDKit (www.rdkit.org) in order to determine category membership. All of this parameters are needed in order to process a chemical of interest through the classification functions. The key/value pairs for each of these categories are as follows:

- DSSTOX_SID: 'dsstox_sid'
- Log Kow: 'logp'
- Molecular Weight: 'mol_weight'
- RDKit.Mol: 'mol'
- SMILES: 'smiles'
- Water Solubility: 'ws', *If using predictions from the expert model suite OPERA (https://ntp.niehs.nih.gov/whatwestudy/niceatm/comptox/ct-opera/opera), these are provided in units as mol/L but the category definitions as implemented rely on units of mg/L.

Chemicals can be provided as a pandas dataFrame containing these columns, a python dictionary containing each of these keys, or a list of dictionaries, each with all of the required keys. In addition, there are built-in error and input checks to warn users if the input type is incompatible with the function being used. 

### Input examples:

```python
from rdkit import Chem
import pandas as pd

# Single-chemical Dictionary
test_chem = {'dsstox_sid': 'DTXSID7020009',
                'smiles': 'CC#N',
                'logp': -0.33971,
                'ws': 12.6417,
                'mol_weight': 41.053,
                'mol': Chem.MolFromSmiles('CC#N')}

# Multi-chemical list of Dictionaries
new_test_chems = [{'dsstox_sid': 'DTXSID3060164',
  'smiles': 'C1=CC=CC=C1C(C1C=CC=CC=1)C1C=CC=CC=1',
  'logp': 5.76,
  'ws': 4.07380277804113e-07,
  'mol_weight': 244.125200512,
  'mol': Chem.MolFromSmiles('C1=CC=CC=C1C(C1C=CC=CC=1)C1C=CC=CC=1')},
 {'dsstox_sid': 'DTXSID7060837',
  'smiles': 'ICCCI',
  'logp': 3.02,
  'ws': 0.0007413102413009177,
  'mol_weight': 295.855896192,
  'mol': Chem.MolFromSmiles('ICCCI')},
 {'dsstox_sid': 'DTXSID9025879',
  'smiles': 'OC(=O)C=CC1C=CC(C=CC(O)=O)=CC=1',
  'logp': 1.99,
  'ws': 0.009120108393559097,
  'mol_weight': 218.0579088,
  'mol': Chem.MolFromSmiles('OC(=O)C=CC1C=CC(C=CC(O)=O)=CC=1')}]

# Multi-chemical Dictionary
test_chems_together = {'dsstox_sid': ['DTXSID3060164','DTXSID7060837'],
  'smiles': ['C1=CC=CC=C1C(C1C=CC=CC=1)C1C=CC=CC=1','ICCCI'],
  'logp': [5.76,3.02],
  'ws': [4.07380277804113e-07,0.0007413102413009177],
  'mol_weight': [244.125200512,295.855896192],
  'mol': [Chem.MolFromSmiles('C1=CC=CC=C1C(C1C=CC=CC=1)C1C=CC=CC=1'),Chem.MolFromSmiles('ICCCI')]}

#Note that the 'mol' attribute can be easily calculated using RDKit with the 'smiles' attribute, so initial data imports likely will not have the 'mol' attribute and instead should be modified as necessary to prep for function use. This process is shown in the final input example below.

# DataFrame
test_chems_df = pd.read_csv('readme_examples.csv', index_col = [0])
test_chems_df['mols'] = [Chem.MolFromSmiles(smile) for smile in test_chems_df['smiles']]

 
```

## User-Focused Functions
epa_ncc include codes that parses the xml file from the OECD Toolbox in order to create many of the tests for the EPA categories. To this end, there are many functions defined that are not meant for general use, rather these permit parsing of the xml and organisation of the various tests and queries. Herein only the functions that allow users to profile chemicals into their respective EPA NCC and to obtain information about those categories are described.

### Code Structure Information
ncc_categories.py defines a Query class. Instances of the Query class correspond to the different available categories. These instances are stored by category key in the dictionary all_tests. Each instance has a .query attribute, which can be applied to an individual chemical in order to obtain a boolean value for whether the given chemical is a member of the specified category. Most queries were built directly from the parsed XML, but some categories required hard-coding of the query tests due to corrupted SMARTS in the XML. These hard-coded categories are italicised in the category list below. 

### Function Definitions
- [**singleQuery**]L1037: A quick method for determining whether a chemical belongs in a specific category.

    - Inputs: 
      - *one_chem*, individual Chemical, provided as a dictionary or DataFrame slice with the keys/columns speficied above. 
      - *category_title*, String representing a category title. Possibilities listed below.
    - Output: *boolean*, value specifies whether x is in Category Title or not

- [**printTree**]L1116: Allows the user to view the testing process for determining whether a chemical belongs in a specific category. Can be run with or without a chemical input.

    - Inputs: 
      - *one_chem*, Default value of x is None but an individual chemical can also be supplied with the same
    constraints as in singleQuery 
      - *category_title*, String representing a category title. Possibilities listed below.
      - *printer*, Boolean with default value True. If True, then this result will be output to the console as
      a print statement. If False, nothing will be printed and the result will instead be a string variable.
    - Output: *printed logic tree*, Each line of the logic tree will contain the query type and all necessary parameters. If data is provided for x, the last value of each line will contain the boolean value for whether x fulfills that piece of the query. 
        - For the XML-originating queries, the first value will be the query ID identifying the query in the XML document. 
        - For hard-coded queries, the first value will instead say CustomQuery and all lines after the first will terminate with "does not process", since the functions for all subqueries are contained within the top branch of the tree only.

If the entire set of logical trees for all chemical categories is desired, add `from epa_ncc.ncc_categories.ncc_categories import all_tests` to your imports and run the following code:

```{python}
the_test_dictionary = {}
for key in all_tests.keys():
    the_test_dictionary[key] = printTree(key,None, printer = False)
```
Now, the variable the_test_dictionary stores strings representing the tests required for classifying chemicals into each category, stored by category title.

- [**queryAll**]L1034: Given a set of chemical(s), returns a DataFrame containing one column for chemical DSSTOXSIDs and individual columns for every category included in all_tests. These columns will contain boolean values, thus describing category membership for the chemical set in a fingerprint-like way. 

    - Inputs: 
      - *chemicals*, A DataFrame, Dictionary, or list of Dictionaries of Chemicals and their attributes, including dsstox_sid, smiles, logp, ws, mol_weight, and RDKIT MolfromSmiles (labelled as 'mol'). There must be keys or column names to match each of these attribute titles.
      - *boolean_outputs*, Default value is False. This function will, by default, output category_df with binary values descripbing category membership. If desired, this matrix can instead be output with boolean values by setting boolean_outputs to True. 
    
    - Output: *category_df*, A DataFrame of chemicals and their category memberships, with an example depicted below:

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>chemicals</th>
      <th>Acid Chlorides</th>
      <th>Acrylamides</th>
      <th>Acrylates/Methacrylates (Acute toxicity)</th>
      <th>Aldehydes (Acute toxicity)</th>
      <th>Aliphatic Amines</th>
      <th>Aluminum Compounds</th>
      <th>Anilines (Acute toxicity)</th>
      <th>Azides (Acute toxicity)</th>
      <th>Benzotriazoles (Acute toxicity)</th>
      <th>...</th>
      <th>Organotins (Chronic toxicity)</th>
      <th>Phenols (Chronic toxicity)</th>
      <th>Phosphinate Esters (Chronic toxicity)</th>
      <th>Polynitroaromatics (Chronic toxicity)</th>
      <th>Substituted Triazines (Chronic toxicity)</th>
      <th>Thiols (Chronic toxicity)</th>
      <th>Vinyl Esters (Chronic toxicity)</th>
      <th>Diazoniums (Chronic toxicity)</th>
      <th>Ethylene Glycol Ethers</th>
      <th>Benzotriazoles</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>DTXSID90480751</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>...</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>1</th>
      <td>DTXSID50939730</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>...</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>DTXSID2036405</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>...</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>DTXSID1024835</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>...</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>DTXSID30878870</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>True</td>
      <td>False</td>
      <td>...</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
<p>5 rows × 67 columns</p>
</div>

- [**listCategories**](https://github.com/laxleary/EPA_Categories/blob/689642cd346f8ae27deaa3df7723742bb4083f3d/categories.py#L1024): Given an individual chemical, this function outputs a list of all categories to which the chemical belongs. 

    - Input: *one_chem*, A DataFrame or Dictionary representing a single chemical and its attributes, including dsstox_sid, smiles, logp, ws, mol_weight, and RDKIT MolfromSmiles (labelled as 'mol'). There must be keys or column names to match each of these attribute titles. 
    - Output: *all_cats*, A list of all categories to which the chemical belongs according to the included tests. 


### Included Categories
- Acid Chlorides
- *Acrylamides*
- *Acrylates/Methacrylates (Acute toxicity)*
- *Acrylates/Methacrylates (Chronic toxicity)*
- *Aldehydes (Acute toxicity)*
- *Aldehydes (Chronic toxicity)*
- *Aliphatic Amines*
- *Alkoxysilanes*
- Aluminum Compounds
- *Aminobenzothiazole Azo Dyes*
- Anhydrides, Carboxylic acid
- Anilines (Acute toxicity)
- Anilines (Chronic toxicity)
- *Anionic Surfactants*
- Azides (Acute toxicity)
- Azides (Chronic toxicity)
- Benzotriazole-hindered phenols
- Benzotriazoles (Acute toxicity)
- Benzotriazoles (Chronic toxicity)
- Boron Compounds
- Cationic (quaternary ammonium) surfactants
- Cobalt
- *Dianilines*
- Diazoniums (Acute toxicity)
- Diazoniums (Chronic toxicity)
- Dichlorobenzidine-based Pigments
- Diisocyanates
- *Dithiocarbamates (Acute toxicity)*
- *Dithiocarbamates (Chronic toxicity)*
- *Epoxides*
- Esters (Acute toxicity)
- Esters (Chronic toxicity)
- *Ethylene Glycol Ethers*
- Hindered Amines
- *Hydrazines and Related Compounds*
- *Imides (Acute toxicity)*
- *Imides (Chronic toxicity)*
- Lanthanides or Rare Earth Metals
- *Neutral Organics*
- Nickel Compounds
- *Nonionic Surfactants*
- *Organotins (Acute toxicity)*
- *Organotins (Chronic toxicity)*
- Peroxides
- Phenolphthaleins
- Phenols (Acute toxicity)
- Phenols (Chronic toxicity)
- Phosphates, Inorganic
- Phosphinate Esters (Acute toxicity)
- Phosphinate Esters (Chronic toxicity)
- *Polynitroaromatics (Acute toxicity)*
- *Polynitroaromatics (Chronic toxicity)*
- Rosin
- Soluble complexes of Zinc
- Stilbene, derivatives of 4,4-bis(triazin-2-ylamino)-
- *Substituted Triazines (Acute toxicity)*
- *Substituted Triazines (Chronic toxicity)*
- *Thiols (Acute toxicity)*
- *Thiols (Chronic toxicity)*
- *Triarylmethane Pigments/Dyes with Non-solubilizing Groups*
- Vinyl Esters (Acute toxicity)
- Vinyl Esters (Chronic toxicity)
- Vinyl Sulfones
- Zirconium Compounds
- *beta-Naphthylamines, Sulfonated*



