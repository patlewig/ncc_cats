# Import all dependencies
import pandas as pd
import numpy as np
import sys
import os
from datetime import datetime
from rdkit import Chem
from pathlib import Path

import os

#Set the working directory for the package
data_path = os.path.join(os.path.dirname(__file__), "data/epa_categories.xml")


# Open the XML file in read mode
with open(data_path,'r', encoding = 'utf8') as f:
    xml=f.read()

# Collapse all XML to a single line (for ease of reading?)
xml=xml.replace('\n','')

# Parse the XML as a tree
import xml.etree.ElementTree as ET
e=ET.parse(data_path).getroot()

# Map children to parents in a dictionary, so that upper elements of the tree can 
# be called by lower elements. This is a many to one mapping. 
parent_map = {c:p for p in e.iter() for c in p}

# Dictionaries converting operator words to mathematical operators and standardizing
# the wording of property names
import operator as op
op_dict={
    'GreaterThan': op.gt,
    'GreaterThanOrEqualTo': op.ge,
    'LessThan': op.lt,
    'LessThanOrEqualTo': op.le
}
prop_dict={
    'log Kow':'logp',
    'Molecular Weight':'mol_weight',
    'Molecular weight':'mol_weight',
    'Water Solubility': 'ws'
}


def define_smart_match(smart):
    """A function that takes in a SMART pattern and outputs a function that identifies whether an input chemical contains that pattern
    Args:
        smart (str): SMARTS pattern as str
    Returns:
        smart_match: Boolean result if SMARTS is a match
    """
    
    pattern=Chem.MolFromSmarts(smart)
    if not pattern:
        return None
    def smart_match(x):
        mol=x['mol']
        ret=True if mol.GetSubstructMatches(pattern) else False
        return ret
    return smart_match


def define_compare(prop,operand,value):
    """A function that takes in a property, operand, and comparative value and outputs a function to determine whether the given property for an input chemical obeys the inequality or equality"""
    def compare(x):
        ret = op_dict[operand](x[prop_dict[prop]],value)
        return ret
    return compare

# I moved these because I think dependencies should be called before they appear to take advantage
# of the Python interpreter
import re
import ast


# Defining the Query class, on which the following functions can be run
class Query:
    
    # Objects in the Query class will have xml, id, logic, subqueries, and category properties, but
    # only the xml (Any Type) is required to create an instance. Some of these will be reset by
    # calling class functions on the instance. 
    def __init__(self,xml,qid=None):
        """ Inputs:
        - xml: The object on which to build an instance of the Class 
        - qid: The numerical ID assigned to the given query, will be set later to match source XML
            \n
        Properties: \n
        - logic: The logic relating the query to other queries
        - subqueries: Queries required to satisfy the query
        - category: The chemical category for which this instance is a query"""
        self.xml=xml
        self.id=qid
        self.logic=None
        self.subqueries=[]


    def write_query(self,qtype:str,tree:dict ):
        """This function creates queries from XML for a Class instance, with inputs of a Query Type (qtype) and 
        a tree. The tree input will be used for logical relations between queries to establish how 
        queries and subqueries relate."""
        # Set the instance type to qtype
        self.type=qtype

        # Check whether the query type is Structural, a Parameter/Property query, or a Logical query
        if qtype=='b:StructureQuery':
            # Pull out the text for the query
            qstring=self.xml.find('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Queries}ComplexSearch').text
            
            # Replace boolean strings with their Python-compatible versions
            qstring=re.sub('false','False',qstring)
            qstring=re.sub('true','True',qstring)
            
            # Interpret the string as Python, save a dictionary of the query pieces
            qdict=ast.literal_eval(qstring)

            # Save the SMART as a variable and set this for the class instance
            smart=qdict['queries'][0]['smart']
            self.smart=smart

            # [Ch3,#1] appears to imply multiple pieces of the SMART that can be matched, so 
            if '[Ch3,#1]' in self.smart:
                split=re.search(r'(.*)\[([^\(\)]*),([^\(\)].*)\]$',self.smart)
                split1=split.group(1)
                split2=split.group(1)+'['+split.group(2)+']'
                smart_match1=define_smart_match(split1)
                smart_match2=define_smart_match(split2)
                def smart_match(x):
                    """Given an input chemical, does the chemical contain either of SMART split 1 or 2?"""
                    return any([smart_match1(x),smart_match2(x)])
            # Without the Ch3, #1 piece just make a function to match the SMART
            else:
                smart_match=define_smart_match(smart)

            # Set the query for the instance to matching the SMART 
            self.query=smart_match


        
        elif qtype=='b:ParameterQuery':
            # Find the operand, property, and value for define_compare
            self.operand=self.xml.find('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Queries}Operand').text
            self.prop=self.xml.find('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Queries}ParameterName').text
            self.value=float(self.xml.find('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Queries}Value').text)
            
            #Create the function needed to apply this query to a chemical
            compare=define_compare(self.prop,self.operand,self.value)

            #Set the query for the instance to checking the property constraints
            self.query=compare


        elif qtype=='LogicalQuery':

            # These capture the logical connectives between other queries, such as "this property AND that property"
            self.logic=self.xml.find('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Engine}Logic').text
            elements=self.xml.find('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Engine}Elements')
            node_ids=[elem.attrib['{http://schemas.microsoft.com/2003/10/Serialization/}Ref']\
                      for elem in elements.findall('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Engine}Query')\
                      if '{http://schemas.microsoft.com/2003/10/Serialization/}Ref' in elem.attrib]
            # Set the query for the instance to checking that the desired property is NOT true
            if self.logic=='Not':
                node_id=node_ids[0] #Should only be one
                sq=tree[node_id]
                self.subqueries=[sq]
                def func(x):
                    return not(sq.query(x))
                self.query=func
            # Set the query for the instance to checking that all desired properties are true
            elif self.logic=='And':
                sqs=[tree[node_id] for node_id in node_ids]
                self.subqueries=sqs
                def func(x):
                    return all([sq.query(x) for sq in self.subqueries])
                self.query=func
            #Set the query for the instance to checking that some property in an OR list is true
            else:
                sqs=[tree[node_id] for node_id in node_ids]
                self.subqueries=sqs
                for orquery in elements.findall('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Engine}Query'):
                    if '{http://www.w3.org/2001/XMLSchema-instance}type' in orquery.attrib:
                        extra_sq=Query(orquery)
                        extra_sq.write_query('b:StructureQuery',tree)
                        sqs.append(extra_sq)      
                def func(x):
                    return any([sq.query(x) for sq in self.subqueries])
                self.query=func
    
    def print_tree(self,x,tabs=0, printer = True):
        """ Given a Class instance and chemical, output the results of applying each query on the chemical
        to the console. For queries with subqueries, the subqueries will be displayed below the query in 
        indented lists. Can also be used to view query conditions without a chemical input."""
        qinfo=(self.id,self.type)
        if self.type=='b:StructureQuery':
            qinfo=qinfo+(self.smart,)
        elif self.type=='b:ParameterQuery':
            qinfo=qinfo+(self.prop,self.value,self.operand)
        elif self.type=='LogicalQuery':
            qinfo=qinfo+(self.logic,)
        try:
            qinfo=qinfo+(bool(self.query(x)),)
        except:
            qinfo=qinfo+('does not process',)
        master_string = '\t'*tabs+str(qinfo)
        if printer:
            print('\t'*tabs+str(qinfo))
        else:
            return master_string
        for sq in self.subqueries:
            sq.print_tree(x,tabs+1, printer)
        if not printer:
            master_string += '\n'+ sq.print_tree(x,tabs+1, printer)
        
    
    

all_tests={}
print_tests = {}
bad_smarts=set()
bad_cats=set()

# This iterates through the XML until it finds a chemical category. Then, for each category, it identifies
#all relevant queries and the contents of those queries.
for elem in e.iter('{http://schemas.microsoft.com/2003/10/Serialization/Arrays}anyType'):
    category=elem.find('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Engine}Caption').text
    queries=elem.find('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Engine}Expression')\
        .find('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Engine}Queries')\
        .findall('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Engine}Query')
    contents=[query.find('{http://schemas.datacontract.org/2004/07/LMC.Profiling.Engine}Content') for query in queries]
    #print(query)
    query_tree={}

    # For each query required by a category... 
    for query in contents:
        attributes=query.attrib
        if '{http://schemas.microsoft.com/2003/10/Serialization/}Id' not in attributes:
            continue

        # Find the qid
        query_id=attributes['{http://schemas.microsoft.com/2003/10/Serialization/}Id']

        # Find the query_type
        query_type=attributes['{http://www.w3.org/2001/XMLSchema-instance}type']

        # Establish a Query instance for the query using its qid
        q=Query(query,query_id)

        # Set the category according to the chemical category this query comes from
        q.category=category

        # Create the query function for this query and print the entire tree related to the query (if it is a 
        # related to a logical query)
        try:
            q.write_query(query_type,query_tree)
        except:
            pass
        # q.query should be a function defining the query if nothing failed, thus if this variable is
        # empty, the category has a problematic query in it. Output the category as a problematic category
        # needing further attention 
        if not q.query or not all([sq.query for sq in q.subqueries]): #Smarts did not compile, sqs needed bc of hidden sqs in or queries
            bad_cats.add(category)

            # In addition, if the query causing the issue was structural, something went wrong in the SMARTS.
            # Store the SMARTs as one that should be looked at. 
            if q.type=='b:StructureQuery':
                bad_smarts.add(q.smart)
        # Store the query in the query tree so that the tree can be built for the entire category   
        query_tree[query_id]=q
    # Store all queries for this category as a tree. The last query called will be the topmost level
    # of the tree, so only that id needs to be preserved     
    all_tests[category]=query_tree[query_id]
    # Store the query id and query type in a dictionary for quick test viewing later
    print_tests[category] = [query_id, query_type]#Final one should always be the top level query hopefully

#The next section contains the hard-coded tests fixed by George

def humanBuiltQuery(function = None, query_words = {'qtype':None, 'smart':None, 'prop':None, 'operand':None, 'logic': None, 'subqueries':None, 'set':None}):
    """This function is an attempt to automate building queries for the hard-coded tests below such that
    print_tree functionality should still work with these newer test types. The query_words dictionary MUST contain
    the key 'qtype' with value 'b:StructureQuery', 'b:ParameterQuery', 'b:ExclusionQuery', or 'LogicalQuery'. Other
    keys should contain the necessary information for print_tree for the given qtype. The function will be the
    query function, as built below. 
    
    Required keys by qtype:
    - Structure or Exclusion Query - 'smart'
    - Parameter Query - 'prop', 'operand', 'value'
    - Logical Query - 'logic', 'subqueries': provided as a list of dictionaries obeying input rules for this function
    """
    new_query = Query(query_words)
    new_query.query = function
    new_query.id = 'CustomQuery'
    qtype = query_words['qtype']
    new_query.type = qtype
    if qtype in ['b:StructureQuery', 'b:ExclusionQuery']:
        new_query.smart = query_words['smart']
    elif qtype == 'b:ParameterQuery':
        new_query.prop = query_words['prop']
        new_query.operand = query_words['operand']
        new_query.value = query_words['value']
    elif qtype == 'LogicalQuery':
        new_query.logic = query_words['logic']
        new_query.subqueries = [humanBuiltQuery(None, query_words['subqueries'][i]) for i in range(len(query_words['subqueries']))]

    return new_query

#Aliphatic amines
def create_test():
    primamine=Chem.MolFromSmarts('[NX3;H2;!$(NC=[O,N,S]);!$(NCN)][CX3]')
    secamine=Chem.MolFromSmarts('[NX3;H1;!$(NC=[O,N,S]);!$(NCN)](C)[CX3]')
    tertamine=Chem.MolFromSmarts('[N;!$(NC=[O,N,S]);!$(NCN)](C)(C)[CX3]')
    def test(x):
        mol=x['mol']
        smiles=x['smiles']
        mw=x['mol_weight']
        return 'c' not in smiles and mw<1000 and '1' not in smiles and (mol.HasSubstructMatch(primamine) or mol.HasSubstructMatch(secamine)\
        or mol.HasSubstructMatch(tertamine)) 
    return test
aa_words = {'qtype':'LogicalQuery', 'logic':'And', \
             'subqueries': [{'qtype':'b:ExclusionQuery', 'smart':'c'}, \
                            {'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value': 1000}, \
                            {'qtype': 'LogicalQuery', 'logic':'or', \
                              'subqueries': [{'qtype': 'b:StructureQuery', 'smart':'[NX3;H2;!$(NC=[O,N,S]);!$(NCN)][CX3]'}, \
                                             {'qtype': 'b:StructureQuery', 'smart':'[NX3;H1;!$(NC=[O,N,S]);!$(NCN)](C)[CX3]'}, \
                                             {'qtype': 'b:StructureQuery', 'smart':'[N;!$(NC=[O,N,S]);!$(NCN)](C)(C)[CX3]'}]}]} 
all_tests['Aliphatic Amines']=humanBuiltQuery(create_test(), aa_words)

#Alkoxysilanes
def create_test():
    alkoxy=Chem.MolFromSmarts('[CX4]O[SiX4]')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        return mw<1000 and mol.HasSubstructMatch(alkoxy)
    return test
alk_words = {'qtype':'LogicalQuery', 'logic':'And', \
             'subqueries': [{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000}, \
                            {'qtype': 'StructureQuery', 'smart':'[CX4]O[SiX4]'}]}
all_tests['Alkoxysilanes']=humanBuiltQuery(create_test(), alk_words)

#Aminobenzothiazole Azo Dyes
def create_test():
    azodye=Chem.MolFromSmiles('N=NC1=NC2=C(S1)C=CC=C2')
    def test(x):
        mol=x['mol'] 
        return mol.HasSubstructMatch(azodye)
    return test
ami_words = {'qtype': 'b:StructureQuery', 'smart':'N=NC1=NC2=C(S1)C=CC=C2'}
all_tests['Aminobenzothiazole Azo Dyes']=humanBuiltQuery(create_test(), ami_words)

#Anionic Surfactants
def has_branching_with_carbon(mol):
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "C" and atom.GetDegree() > 2:
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            # print(f"Carbon index {atom.GetIdx()} degree: {atom.GetDegree()}, neighbors: {neighbors}")
            # Only count branching caused by carbon neighbors
            carbon_neighbors = [n for n in neighbors if n == "C"]
            if len(carbon_neighbors) > 2:
                # print(f"Branching detected at carbon index {atom.GetIdx()} with neighbors {neighbors}")
                return True
    return False
from rdkit.Chem import rdmolops

def is_straight_alkyl_chain(mol):
    carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() == "C"]
    for idx in carbons:
        atom = mol.GetAtomWithIdx(idx)
        # Ensure no carbon has more than 2 non-functional group bonds
        if atom.GetDegree() > 2:
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if not all(n in {"O", "S"} for n in neighbors if n != "C"):
                return False
    return True

def anionic_test():
    sulfate = Chem.MolFromSmarts('COS(=O)(=O)[OH,O-]')
    sulfonate = Chem.MolFromSmarts('CS(=O)(=O)[OH,O-]')
    phosphate = Chem.MolFromSmarts('COP([OH1])([OH1])=O')
    carboxylic = Chem.MolFromSmarts('[CX3;!$(Cc)](=O)[OX2H1]')
    silicic = Chem.MolFromSmarts('[Si][OX2H]')
    
    if not all([sulfate, sulfonate, phosphate, carboxylic, silicic]):
        raise ValueError("One or more SMARTS patterns failed to initialize.")
    
    def test(x):
        mol = x['mol']
        
        # Character filtering
        if mol is None:
            return False
        
        # Exclude branching with carbon atoms
        
        if has_branching_with_carbon(mol):
            return False
        
        # Ensure a straight alkyl chain
        
        if not is_straight_alkyl_chain(mol):
            return False
        
        # Substructure matching
        return (
            mol.HasSubstructMatch(sulfate) or 
            mol.HasSubstructMatch(sulfonate) or 
            mol.HasSubstructMatch(phosphate) or 
            mol.HasSubstructMatch(carboxylic) or 
            mol.HasSubstructMatch(silicic)
        )
    
    return test


ani_words = {'qtype': 'LogicalQuery', 'logic':'And', \
             'subqueries': [{'qtype':'b:StructureQuery', 'smart':'Straight Alkyl Chain - Cs are adjacent in SMILES'}, \
                            {'qtype':'b:StructureQuery', 'smart':'Has branching with Carbon'},\
                            {'qtype':'LogicalQuery', 'logic':'Or', \
                             'subqueries': [{'qtype':'b:StructureQuery', 'smart':'COS(=O)(=O)[OH,O-]'},\
                                            {'qtype':'b:StructureQuery', 'smart':'CS(=O)(=O)[OH,O-]'}, \
                                            {'qtype':'b:StructureQuery', 'smart':'COP([OH1])([OH1])=O'}, \
                                            {'qtype':'b:StructureQuery', 'smart':'[CX3;!$(Cc)](=O)[OX2H1]'}, \
                                            {'qtype':'b:StructureQuery', 'smart':'[Si][OX2H]'}]}]}
all_tests['Anionic Surfactants']=humanBuiltQuery(anionic_test(), ani_words)

# #Benzotriazoles
# Original code has this test twice, but this one would be overwritten by the later copy
# so it is commented out for now to avoid confusion when searching through the hard-coded tests
# def create_test():
#     benzotriazole=Chem.MolFromSmarts('n1c2ccccc2nn1')
#     def test(x):
#         mol=x['mol']
#         return mol.HasSubstructMatch(benzotriazole)
#     return test
# benzo_words = {'qtype': 'b:StructureQuery', 'smart':'n1c2ccccc2nn1'}
# all_tests['Benzotriazoles']=humanBuiltQuery(create_test(), benzo_words)

#Dianilines
def create_test():
    dianiline=Chem.MolFromSmarts('c1cc([NH2])ccc1[CH2,O,N,S]c1ccccc1')
    not_dianiline1=Chem.MolFromSmarts('c1ccccc1[A]~[A]')
    not_dianiline2=Chem.MolFromSmarts('c1ccccc1[A](c)c')
    def test(x):
        mol=x['mol']
        return not mol.HasSubstructMatch(not_dianiline1) and not mol.HasSubstructMatch(not_dianiline2)\
        and len(mol.GetSubstructMatches(dianiline))==2 #lol
    return test
dia_words = {'qtype': 'LogicalQuery', 'logic':'And', \
              'subqueries': [{'qtype':'b:StructureQuery', 'smart':'c1cc([NH2])ccc1[CH2,O,N,S]c1ccccc1, Needs Two Copies'},\
                             {'qtype':'b:ExclusionQuery', 'smart':'c1ccccc1[A]~[A]'},\
                             {'qtype':'b:ExclusionQuery', 'smart':'c1ccccc1[A](c)c'}]}
all_tests['Dianilines']=humanBuiltQuery(create_test(), dia_words)

#Dithiocarbonates 
def create_test():
    ethylenebisdithiocarbamate=Chem.MolFromSmiles('SC(=S)NCCNC(=S)S')
    dithiocarbamates=[]
    for i in range(1,5):
        for j in range(1,5):
            mol=Chem.MolFromSmiles('C'*i + 'NC(=S)S' + 'C'*j)
            dithiocarbamates.append(mol)
    dithiocarbamates.append(ethylenebisdithiocarbamate)
    def test(x):
        mol=x['mol']
        return x['mol_weight']<1000 and x['logp']<5 and any([mol.HasSubstructMatch(dithiocarbamate) for dithiocarbamate in dithiocarbamates])
    return test
dithiocarbamate_list = ['CNC(=S)SC',
                    'CCSC(=S)NC',
                    'CCCSC(=S)NC',
                    'CCCCSC(=S)NC',
                    'CCNC(=S)SC',
                    'CCNC(=S)SCC',
                    'CCCSC(=S)NCC',
                    'CCCCSC(=S)NCC',
                    'CCCNC(=S)SC',
                    'CCCNC(=S)SCC',
                    'CCCNC(=S)SCCC',
                    'CCCCSC(=S)NCCC',
                    'CCCCNC(=S)SC',
                    'CCCCNC(=S)SCC',
                    'CCCCNC(=S)SCCC',
                    'CCCCNC(=S)SCCCC',
                    'S=C(S)NCCNC(=S)S']
dith_words1 = {'qtype':'LogicalQuery', 'logic':'And', \
              'subqueries':[{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},\
                            {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThan', 'value':5}, \
                            {'qtype':'LogicalQuery', 'logic':'Or', 'subqueries':[]}]}
for smart in dithiocarbamate_list:
    dith_words1['subqueries'][2]['subqueries'].append({'qtype':'b:StructureQuery', 'smart':smart})
all_tests['Dithiocarbamates (Acute toxicity)']=humanBuiltQuery(create_test(), dith_words1)
def create_test():
    ethylenebisdithiocarbamate=Chem.MolFromSmiles('SC(=S)NCCNC(=S)S')
    dithiocarbamates=[]
    for i in range(1,5):
        for j in range(1,5):
            mol=Chem.MolFromSmiles('C'*i + 'NC(=S)S' + 'C'*j)
            dithiocarbamates.append(mol)
    dithiocarbamates.append(ethylenebisdithiocarbamate)
    def test(x):
        mol=x['mol']
        return x['mol_weight']<1000 and x['logp']>=5 and x['logp']<19 and any([mol.HasSubstructMatch(dithiocarbamate) for dithiocarbamate in dithiocarbamates])
    return test
dith_words2 = {'qtype':'LogicalQuery', 'logic':'And', \
              'subqueries':[{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},\
                            {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'GreaterThanOrEqualTo', 'value':5},\
                            {'qtype':'b:ParameterQuery', 'prop': 'Log Kow', 'operand':'LessThan', 'value':19}, \
                            {'qtype':'LogicalQuery', 'logic':'Or', 'subqueries':[]}]}
for smart in dithiocarbamate_list:
    dith_words2['subqueries'][3]['subqueries'].append({'qtype':'b:StructureQuery', 'smart':smart})
all_tests['Dithiocarbamates (Chronic toxicity)']=humanBuiltQuery(create_test(), dith_words2)

#Ethylene Glycol Ethers
#Have to enumerate       
def create_test():
    #match_mols and phenyl_mols do not appear to get used in the actual test?
    match_mols=[]
    for i in range(1,8):
        for j in range(0,8):
            for k in range(1,4):
                smart='C'*i+'OCC'*k+'O'+'C'*j
                match_mols.append(Chem.MolFromSmiles(smart))
    phenyl_mols=[]
    for i in range(0,7):
        for k in range(1,3):
            for l in range(0,3): #Technically could be any number but this is difficult to implement
                phenyl_smart='c1ccccc1'+'C'*l+'OCC'*k+'O'+'C'*i
                phenyl_mols.append(Chem.MolFromSmiles(phenyl_smart))
    
    def test(x):
        smiles=x['smiles']
        if set(smiles)-set(['C','c','O','1','(',')']):
            return False
        if smiles.count('O')<2:
            return False
        os=[i for i,o in enumerate(smiles) if o=='O']
        between_os=[smiles[(start+1):end] for start,end in zip(os,os[1:])]
        if any([between!='CC' for between in between_os]):
            return False
        m=re.compile('1.*O.*1')
        if m.findall(smiles):
            return False
        carbon1=smiles[0:min(os)]
        carbon2=smiles[(max(os)+1):]
        if carbon1.count('C')>7 or carbon2.count('C')>7:
            return False
        if carbon1.count('c')>6 or carbon2.count('c')>6:
            return False
        if not carbon1 and not carbon2:
            return False
        else:
            return True
    return test
ethyl_words = {'qtype':'LogicalQuery', 'logic':'And',
               'subqueries':[{'qtype':'b:StructureQuery', 'smart':"All symbols contained in ['C','c','O','1','(',')']"}, 
                             {'qtype':'b:StructureQuery', 'smart':'At Least Two O'},
                             {'qtype':'b:StructureQuery', 'smart':'No two Cs between Os'},
                             {'qtype':'b:StructureQuery', 'smart':'No repeating chains bewteen 1-O-1'},
                             {'qtype':'b:StructureQuery', 'smart':'There are fewer than 7 Cs before the first O and after the last O'},
                             {'qtype':'b:StructureQuery', 'smart':'There are fewer than 7 cs before the first O and after the last O'},
                             {'qtype':'b:StructureQuery', 'smart':'There is something before the first O or after the last O'}]}
all_tests['Ethylene Glycol Ethers']=humanBuiltQuery(create_test(), ethyl_words)

#Neutral Organics
#Verhaar scheme, see paper called Classifying Environmental Pollutants
def create_test():
    def test(x):
        mol=x['mol']
        if mol.HasSubstructMatch(Chem.MolFromSmarts('[!C;!c;!N;!O;!F;!Cl;!Br,I]')): #Rule 0.1 and 1.1
            return False
        logp=x['logp']
        if logp>8:
            return False
        mw=x['mol_weight']
        if mw>1000:
            return False
        if not mol.HasSubstructMatch(Chem.MolFromSmarts('[!C;!c]')): #Rule 1.3
            return True
        elif not mol.HasSubstructMatch(Chem.MolFromSmarts('[!C;!c;!Cl;!Br;!F]'))\
        and not mol.HasSubstructMatch(Chem.MolFromSmarts('[Cl,Br,F]C[$(C=C),$(Cc)]')): #Rule 1.4
            return True
        elif not mol.HasSubstructMatch(Chem.MolFromSmarts('[!C;!c;!O;!Cl;!Br;!F]')): #Rule 1.5
            if mol.HasSubstructMatch(Chem.MolFromSmarts('COC'))\
            and not mol.HasSubstructMatch(Chem.MolFromSmarts('COOC'))\
            and not mol.HasSubstructMatch(Chem.MolFromSmarts('C1OC1')): #Rule 1.5.1 and 1.7
                return True
            elif mol.HasSubstructMatch(Chem.MolFromSmarts('[C;!$(C=O)][OH]'))\
            and not mol.HasSubstructMatch(Chem.MolFromSmarts('C=CCO'))\
            and not mol.HasSubstructMatch(Chem.MolFromSmarts('C#CCO'))\
            and not mol.HasSubstructMatch(Chem.MolFromSmarts('cCO')): #Rule 1.5.2, 1.5.3, and 1.7CCCCOCCOCCO
                return True
            elif mol.HasSubstructMatch(Chem.MolFromSmarts('[C;!$(CO)]=O'))\
            and not mol.HasSubstructMatch(Chem.MolFromSmarts('[$(cC),$(C=C)]C=O'))\
            and not mol.HasSubstructMatch(Chem.MolFromSmarts('[Cl,Br]C=O'))\
            and not mol.HasSubstructMatch(Chem.MolFromSmarts('[Cl,Br]CC=O')): #Rule 1.5.4 and 1.7
                return True
            else:
                return False
        elif not mol.HasSubstructMatch(Chem.MolFromSmarts('[!C;!N]'))\
        and mol.HasSubstructMatch(Chem.MolFromSmarts('C[NH,NH0]')): #Rule 1.6
            return True
        else: 
            return False
    return test
neutral_words = {'qtype':'LogicalQuery', 'logic':'This test contains more than 21 branches following the Verhaar Scheme from: https://doi.org/10.1016/0045-6535(92)90280-5', 'subqueries':[]}
all_tests['Neutral Organics']=humanBuiltQuery(create_test(), neutral_words)

#Nonionic Surfactants

#This is immediately repeated and overwritten, so this copy is commented out for clarity on which copy is used

# nonsurf1=Chem.MolFromSmarts('COCCO')
# nonsurf2=Chem.MolFromSmarts('COCCOC')
# def test(x):
#     mol=x['mol']
#     return mol.HasSubstructMatch(nonsurf1) or mol.HasSubstructMatch(nonsurf2)
# import re
# def test(x):
#     smiles=x['smiles']
#     if '(' in smiles:
#         return False
#     split_smiles=smiles.split('O')
#     if len(split_smiles)==1:
#         return False
#     mol=x['mol']
#     if not mol.HasSubstructMatch(Chem.MolFromSmiles('COC')) or mol.HasSubstructMatch(Chem.MolFromSmiles('C=O')):
#         return False
#     return not any([re.search(r'[^C]',c) for c in split_smiles])
# all_tests['Nonionic Surfactants']=test

#Nonionic Surfactants
import math
def create_test():
    def test(x):
        mol=x['mol']
        atoms=[a for a in x['smiles'].lower() if a.isalpha()]
        if set(atoms)-set(['c','o']):
            return False
        return mol.HasSubstructMatch(Chem.MolFromSmarts('[CH3][CR0][CR0][CR0][CR0][CR0]')) and\
        atoms.count('o')>1 and\
        (math.floor(len(mol.GetSubstructMatches(Chem.MolFromSmarts('O[CH2][CH2]')))/2)+1)==len(mol.GetSubstructMatches(Chem.MolFromSmarts('[O]')))
    return test
nonio_words = {'qtype':'LogicalQuery', 'logic':'And',
               'subqueries':[{'qtype':'b:StructureQuery', 'smart':'[CH3][CR0][CR0][CR0][CR0][CR0]'},
                             {'qtype':'b:StructureQuery', 'smart':'There is more than 1 o'},
                             {'qtype':'b:StructureQuery', 'smart':'The number of O[CH2][CH2] divided by 2, plus 1, rounded down, is the number of [O]'}]}
all_tests['Nonionic Surfactants']=humanBuiltQuery(create_test(), nonio_words)

#Organotins (Acute toxicity) and Organotins (Chronic toxicity)

def create_test():
    organotin=Chem.MolFromSmarts('C[Sn]') 
    def test(x):
        mol=x['mol']
        return x['mol_weight']<1000 and mol.HasSubstructMatch(organotin) and x['logp']<=13.7
    return test
orga_words = {'qtype':'LogicalQuery', 'logic':'And',
              'subqueries':[{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},
                            {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThanOrEqualTo', 'value':13.7},
                            {'qtype':'b:StructureQuery', 'smart':'C[Sn]'}]}
all_tests['Organotins (Acute toxicity)']=humanBuiltQuery(create_test(), orga_words)
def create_test():
    organotin=Chem.MolFromSmarts('C[Sn]') 
    def test(x):
        mol=x['mol']
        return x['mol_weight']<1000 and mol.HasSubstructMatch(organotin) and x['logp']>=13.7
    return test
orga_words2 = {'qtype':'LogicalQuery', 'logic':'And',
              'subqueries':[{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},
                            {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'GreaterThanOrEqualTo', 'value':13.7},
                            {'qtype':'b:StructureQuery', 'smart':'C[Sn]'}]}
all_tests['Organotins (Chronic toxicity)']=humanBuiltQuery(create_test(), orga_words2)

#Polynitroaromatics (Acute toxicity) and Polynitroaromatics (Chronic toxicity)
#MW < 1000

def create_test():
    polynitroaromatic=Chem.MolFromSmarts('ON(=O)[$(c1c(N(O)=O)cccc1),$(c1cc(N(O)=O)ccc1),$(c1ccc(N(O)=O)cc1),$(c1cncc(N(O)=O)c1)]')
    def test(x):
        mol=x['mol']
        return x['mol_weight']<1000 and mol.HasSubstructMatch(polynitroaromatic) and x['logp']<7
    return test
polyn_words = {'qtype':'LogicalQuery', 'logic':'And',
               'subqueries':[{'qtype':'b:StructureQuery', 'smart':'ON(=O)[$(c1c(N(O)=O)cccc1),$(c1cc(N(O)=O)ccc1),$(c1ccc(N(O)=O)cc1),$(c1cncc(N(O)=O)c1)]'},
                             {'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},
                             {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThan', 'value':7}]}
all_tests['Polynitroaromatics (Acute toxicity)']=humanBuiltQuery(create_test(), polyn_words)
def create_test():
    polynitroaromatic=Chem.MolFromSmarts('N[$(c1c(N)cccc1),$(c1cc(N)ccc1),$(c1ccc(N)cc1),$(c1cncc(N)c1)]')
    def test(x):
        mol=x['mol']
        return x['mol_weight']<1000 and mol.HasSubstructMatch(polynitroaromatic) and x['logp']>=10
    return test
polyn_words2 = {'qtype':'LogicalQuery', 'logic':'And',
               'subqueries':[{'qtype':'b:StructureQuery', 'smart':'N[$(c1c(N)cccc1),$(c1cc(N)ccc1),$(c1ccc(N)cc1),$(c1cncc(N)c1)]'},
                             {'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},
                             {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'GreaterThanOrEqualTo', 'value':10}]}
all_tests['Polynitroaromatics (Chronic toxicity)']=humanBuiltQuery(create_test(), polyn_words2)

#Substituted Triazines (Acute toxicity) and Substituted Triazines (Chronic toxicity)
#logp<5
#MW<1000
def create_test():
    subtriazine=Chem.MolFromSmarts('[$(n1nnccc1.[!#1]),$(n1ncncc1.[!#1]),$(n1cncnc1.[!#1])]')#[!H] did not work as expected with aromatics
    def test(x):
        mol=x['mol']
        return x['mol_weight']<1000 and mol.HasSubstructMatch(subtriazine) and x['logp']<5
    return test
subtri_words1 = {'qtype':'LogicalQuery', 'logic':'And', \
                 'subqueries':[{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000}, \
                               {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThan', 'value':5},\
                               {'qtype':'b:StructureQuery', 'smart':'[$(n1nnccc1.[!#1]),$(n1ncncc1.[!#1]),$(n1cncnc1.[!#1])]'}]}
all_tests['Substituted Triazines (Acute toxicity)']=humanBuiltQuery(create_test(), subtri_words1)
def create_test():
    subtriazine=Chem.MolFromSmarts('[$(n1nnccc1.[!#1]),$(n1ncncc1.[!#1]),$(n1cncnc1.[!#1])]')#[!H] did not work as expected with aromatics
    def test(x):
        mol=x['mol']
        return x['mol_weight']<1000 and mol.HasSubstructMatch(subtriazine) and x['logp']>5 and x['logp']<=8
    return test
subtri_words2 = {'qtype':'LogicalQuery', 'logic':'And', \
                 'subqueries':[{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000}, \
                               {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'GreaterThan', 'value':5},\
                               {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThanOrEqualTo', 'value':8},\
                               {'qtype':'b:StructureQuery', 'smart':'[$(n1nnccc1.[!#1]),$(n1ncncc1.[!#1]),$(n1cncnc1.[!#1])]'}]}
all_tests['Substituted Triazines (Chronic toxicity)']=humanBuiltQuery(create_test(), subtri_words2)

def convert_ppb(x): #OPERA results stored as mol/L
    ws=x['ws']
    mol_weight=x['mol_weight']
    return ws*mol_weight*10**3 #Based on what the Toolbox has done - seems that the OECD Toolbox xl has implemented WSs in units of mg/L. OPERA predictions are in log10 mol/L so the conversion needs to be solubility in units of mol/L * MW * 1000

#Triarylmethane Pigments/Dyes with Non-solubilizing Groups
def create_test():
    para_permutations='[NH2,O,$([NH1][CH3]),$([NH1][CH2][CH3]),$(N([CH3])[CH3]),$(N([CH3])[CH2][CH3]),$(N([CH2][CH3])[CH2][CH3])]'
    triphenylmethane=Chem.MolFromSmarts('[cH]1[cH]c({})[cH][cH]c1C(c2[cH][cH]c({})[cH][cH]2)=C3[CH]=[CH]C(=[NH,O])[CH]=[CH]3'.format(para_permutations,para_permutations))
    diphenylnaphthylmethane=Chem.MolFromSmarts('[cH]1[cH]c({})[cH][cH]c1C(c2[cH][cH]c({})[cH]3[cH][cH][cH][cH][cH]32)=C3[CH]=[CH]C(=[NH,O])[CH]=[CH]3'.format(para_permutations,para_permutations))
    def test(x):
        mol=x['mol']
        return convert_ppb(x)>1 and (mol.HasSubstructMatch(triphenylmethane) or (mol.HasSubstructMatch(diphenylnaphthylmethane)))
    return test
triar_words = {'qtype':'LogicalQuery', 'logic':'And', \
               'subqueries':[{'qtype':'b:ParameterQuery', 'prop':'Water Solubility', 'operand':'GreaterThan', 'value':1},\
                             {'qtype':'LogicalQuery', 'logic':'Or', \
                              'subqueries':[{'qtype':'b:StructureQuery', 'smart':'[cH]1[cH]c({})[cH][cH]c1C(c2[cH][cH]c({})[cH][cH]2)=C3[CH]=[CH]C(=[NH,O])[CH]=[CH]3 Formatted with Parapermutations: [NH2,O,$([NH1][CH3]),$([NH1][CH2][CH3]),$(N([CH3])[CH3]),$(N([CH3])[CH2][CH3]),$(N([CH2][CH3])[CH2][CH3])]'},\
                                            {'qtype':'b:StructureQuery', 'smart':'[cH]1[cH]c({})[cH][cH]c1C(c2[cH][cH]c({})[cH]3[cH][cH][cH][cH][cH]32)=C3[CH]=[CH]C(=[NH,O])[CH]=[CH]3 Formatted with Parapermutations: [NH2,O,$([NH1][CH3]),$([NH1][CH2][CH3]),$(N([CH3])[CH3]),$(N([CH3])[CH2][CH3]),$(N([CH2][CH3])[CH2][CH3])]'}]}]}
all_tests['Triarylmethane Pigments/Dyes with Non-solubilizing Groups']=humanBuiltQuery(create_test(), triar_words)

#beta-Naphthylamines, Sulfonated
def create_test():
    smarts=[]
    match_mols=[]
    prefix='[NH2]c1[cH1,$(cO)]'
    suffix='[cH][cH]1'
    for c1 in range(1,4):
        for c2 in range(c1+1,5):
            smarts.append(prefix+'c2'+'[cH]'*(c1-1)+'[cH,$(c[OH]),$(c[NH2])]'+'[cH]'*(c2-c1-1)+'c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])'+'[cH]'*(4-c2)+'c2'+suffix)
            smarts.append(prefix+'c2'+'[cH]'*(c1-1)+'c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])'+'[cH]'*(c2-c1-1)+'[cH1,$(c[OH]),$(c[NH2])]'+'[cH]'*(4-c2)+'c2'+suffix)
    match_mols=[Chem.MolFromSmarts(smart) for smart in smarts]
    def test(x):
        mol=x['mol']
        naph_matches=[True for match in match_mols[:] if mol.HasSubstructMatch(match) and match.HasSubstructMatch(mol)]
        return any(naph_matches)
    return test
beta_words = {'qtype':'LogicalQuery', 'logic':'or', \
              'subqueries':[]}
smarts_list = ['[NH2]c1[cH1,$(cO)]c2[cH,$(c[OH]),$(c[NH2])]c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])[cH][cH]c2[cH][cH]1',
 '[NH2]c1[cH1,$(cO)]c2c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])[cH1,$(c[OH]),$(c[NH2])][cH][cH]c2[cH][cH]1',
 '[NH2]c1[cH1,$(cO)]c2[cH,$(c[OH]),$(c[NH2])][cH]c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])[cH]c2[cH][cH]1',
 '[NH2]c1[cH1,$(cO)]c2c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])[cH][cH1,$(c[OH]),$(c[NH2])][cH]c2[cH][cH]1',
 '[NH2]c1[cH1,$(cO)]c2[cH,$(c[OH]),$(c[NH2])][cH][cH]c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])c2[cH][cH]1',
 '[NH2]c1[cH1,$(cO)]c2c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])[cH][cH][cH1,$(c[OH]),$(c[NH2])]c2[cH][cH]1',
 '[NH2]c1[cH1,$(cO)]c2[cH][cH,$(c[OH]),$(c[NH2])]c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])[cH]c2[cH][cH]1',
 '[NH2]c1[cH1,$(cO)]c2[cH]c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])[cH1,$(c[OH]),$(c[NH2])][cH]c2[cH][cH]1',
 '[NH2]c1[cH1,$(cO)]c2[cH][cH,$(c[OH]),$(c[NH2])][cH]c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])c2[cH][cH]1',
 '[NH2]c1[cH1,$(cO)]c2[cH]c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])[cH][cH1,$(c[OH]),$(c[NH2])]c2[cH][cH]1',
 '[NH2]c1[cH1,$(cO)]c2[cH][cH][cH,$(c[OH]),$(c[NH2])]c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])c2[cH][cH]1',
 '[NH2]c1[cH1,$(cO)]c2[cH][cH]c([$(S(=O)(=O)[OH]),$(S(=O)(=O)[CH2][CH2]S[OH3])])[cH1,$(c[OH]),$(c[NH2])]c2[cH][cH]1']
for smart in smarts_list:
    beta_words['subqueries'].append({'qtype':'b:StructureQuery', 'smart':smart + ", Must Match in Both Directions"})
all_tests['beta-Naphthylamines, Sulfonated']=humanBuiltQuery(create_test(), beta_words)

#Aldehydes
#Turns out these are just wrong in the toolbox, although does compile
def create_test():
    formaldehyde=Chem.MolFromSmarts('[CH2](=O)') #Needs to be special case because buggy way RDKit handles hydrogens
    aldehyde=Chem.MolFromSmarts('[CH1](=[O])[C,c]')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        logp=x['logp']
        return (mol.HasSubstructMatch(formaldehyde) or mol.HasSubstructMatch(aldehyde)) and mw<1000 and logp<=6
    return test
alde_words = {'qtype':'LogicalQuery', 'logic':'And', \
              'subqueries':[{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},\
                            {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand': 'LessThanOrEqualTo', 'value':6},\
                            {'qtype':'LogicalQuery', 'logic':'Or', \
                             'subqueries':[{'qtype':'b:StructureQuery', 'smart':'[CH2](=O)'},\
                                           {'qtype':'b:StructureQuery', 'smart':'[CH1](=[O])[C,c]'}]}]}
all_tests['Aldehydes (Acute toxicity)']=humanBuiltQuery(create_test(), alde_words)
def create_test():
    formaldehyde=Chem.MolFromSmarts('[CH2](=O)') #Needs to be special case because buggy way RDKit handles hydrogens
    aldehyde=Chem.MolFromSmarts('[CH1](=[O])[C,c]')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        logp=x['logp']
        return (mol.HasSubstructMatch(formaldehyde) or mol.HasSubstructMatch(aldehyde)) and mw<1000 and logp>6
    return test
alde_words2 = {'qtype':'LogicalQuery', 'logic':'And', \
              'subqueries':[{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},\
                            {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand': 'GreaterThan', 'value':6},\
                            {'qtype':'LogicalQuery', 'logic':'Or', \
                             'subqueries':[{'qtype':'b:StructureQuery', 'smart':'[CH2](=O)'},\
                                           {'qtype':'b:StructureQuery', 'smart':'[CH1](=[O])[C,c]'}]}]}
all_tests['Aldehydes (Chronic toxicity)']=humanBuiltQuery(create_test(), alde_words2)

# #Benzotriazoles
# #Not a valid smarts from toolbox
# def create_test():
#     benzotriazole=Chem.MolFromSmarts('N1N=NC2=C1C=CC=C2')
#     def test(x):
#         mol=x['mol']
#         return mol.HasSubstructMatch(benzotriazole)
#     return test
# benzo_words = {'qtype': 'b:StructureQuery', 'smart':'N1N=NC2=C1C=CC=C2'}
# all_tests['Benzotriazoles']=humanBuiltQuery(create_test(), benzo_words)

#Imides
#Doesn't work if carbons are part of aromatic
def create_test():
    imide=Chem.MolFromSmarts('C(=O)NC(=O)')
    not_imide=Chem.MolFromSmarts('c1C(=O)NC(=O)ccccc1')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        logp=x['logp']
        return mol.HasSubstructMatch(imide) and not mol.HasSubstructMatch(not_imide) and logp<=5 and mw<1000
    return test
imi_words = {'qtype':'LogicalQuery', 'logic':'And', \
             'subqueries':[{'qtype':'b:StructureQuery', 'smart':'C(=O)NC(=O)'},\
                           {'qtype':'b:ExclusionQuery', 'smart':'c1C(=O)NC(=O)ccccc1'},\
                           {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThanOrEqualTo', 'value':5},\
                           {'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000}]}
all_tests['Imides (Acute toxicity)']=humanBuiltQuery(create_test(), imi_words)
def create_test():
    imide=Chem.MolFromSmarts('C(=O)NC(=O)')
    not_imide=Chem.MolFromSmarts('c1(C(=O)NC(=O))ccccc1')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        logp=x['logp']
        return mol.HasSubstructMatch(imide) and not mol.HasSubstructMatch(not_imide) and logp>5 and logp<8 and mw<1000
    return test
imi_words2 = {'qtype':'LogicalQuery', 'logic':'And', \
             'subqueries':[{'qtype':'b:StructureQuery', 'smart':'C(=O)NC(=O)'},\
                           {'qtype':'b:ExclusionQuery', 'smart':'c1C(=O)NC(=O)ccccc1'},\
                           {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'GreaterThan', 'value':5},\
                           {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThan', 'value':8},\
                           {'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000}]}
all_tests['Imides (Chronic toxicity)']=humanBuiltQuery(create_test(), imi_words2)

#Hydrazines and related compounds
def create_test():
    hydra1=Chem.MolFromSmarts('[NX3][NX3]')
    hydra2=Chem.MolFromSmarts('[CX3]=[NX2][NX2]')
    hydra3=Chem.MolFromSmarts('[CX3](=O)[NX2][NX3]')
    hydra4=Chem.MolFromSmarts('[NX2][CX3](=O)[NX2][NX3]')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        return (mol.HasSubstructMatch(hydra1) or mol.HasSubstructMatch(hydra2)\
               or mol.HasSubstructMatch(hydra3) or mol.HasSubstructMatch(hydra4)) and mw<500
    return test
hydra_words = {'qtype':'LogicalQuery', 'logic':'And', \
               'subqueries':[{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':500}, \
                             {'qtype':'LogicalQuery', 'logic':'Or', \
                              'subqueries':[{'qtype':'b:StructureQuery', 'smart':'[NX3][NX3]'},\
                                            {'qtype':'b:StructureQuery', 'smart':'[CX3]=[NX2][NX2]'},\
                                            {'qtype':'b:StructureQuery', 'smart':'[CX3](=O)[NX2][NX3]'},\
                                            {'qtype':'b:StructureQuery', 'smart':'[NX2][CX3](=O)[NX2][NX3]'}]}]}
all_tests['Hydrazines and Related Compounds']=humanBuiltQuery(create_test(), hydra_words)

#Thiols
def create_test():
    thiol=Chem.MolFromSmarts('[C,c][SX2H]')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        logp=x['logp']
        return mol.HasSubstructMatch(thiol) and mw<1000 and logp<6.5
    return test
thio_words = {'qtype':'LogicalQuery', 'logic':'And', \
              'subqueries': [{'qtype':'b:StructureQuery', 'smart':'[C,c][SX2H]'},\
                             {'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},\
                             {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThan', 'value':6.5}]}
all_tests['Thiols (Acute toxicity)']=humanBuiltQuery(create_test(), thio_words)
def create_test():
    thiol=Chem.MolFromSmarts('[C,c][SX2H]')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        logp=x['logp']
        return mol.HasSubstructMatch(thiol) and mw<1000 and logp>=6.5 and logp<9
    return test
thio_words2 = {'qtype':'LogicalQuery', 'logic':'And', \
              'subqueries': [{'qtype':'b:StructureQuery', 'smart':'[C,c][SX2H]'},\
                             {'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},\
                             {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'GreaterThanOrEqualTo', 'value':6.5},\
                             {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThan', 'value':9}  ]}
all_tests['Thiols (Chronic toxicity)']=humanBuiltQuery(create_test(), thio_words2)

#Acrylamides
def create_test():
    acrylamide1=Chem.MolFromSmarts('[CH2]=[CH1]C(=O)[NH,NH2]')
    acrylamide2=Chem.MolFromSmarts('[CH2]=C([CH3])C(=O)[NH,NH2]')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        logp=x['logp']
        return (mol.HasSubstructMatch(acrylamide1) or mol.HasSubstructMatch(acrylamide2)) and mw<1000 and logp<8
    return test
acry_words = {'qtype':'LogicalQuery', 'logic':'And', \
              'subqueries': [{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},\
                             {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThan', 'value':8}, \
                             {'qtype':'LogicalQuery', 'logic':'Or',\
                              'subqueries':[{'qtype':'b:StructureQuery', 'smart':'[CH2]=[CH1]C(=O)[NH,NH2]'},\
                                            {'qtype':'b:StructureQuery', 'smart':'[CH2]=C([CH3])C(=O)[NH,NH2]'}]}]}
all_tests['Acrylamides']=humanBuiltQuery(create_test(), acry_words)

#Acrylates/Methacrylates
def create_test():
    acrylate=Chem.MolFromSmarts('[CH2]=[CH]C(=O)O')
    methacrylate=Chem.MolFromSmarts('[CH2]=C([CH3])C(=O)O')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        logp=x['logp']
        return (mol.HasSubstructMatch(acrylate) or mol.HasSubstructMatch(methacrylate)) and logp<=5 and mw<1000
    return test
acry_words = {'qtype':'LogicalQuery', 'logic':'And', \
              'subqueries': [{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},\
                             {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThanOrEqualTo', 'value':5}, \
                             {'qtype':'LogicalQuery', 'logic':'Or',\
                              'subqueries':[{'qtype':'b:StructureQuery', 'smart':'[CH2]=[CH]C(=O)O'},\
                                            {'qtype':'b:StructureQuery', 'smart':'[CH2]=C([CH3])C(=O)O'}]}]}
all_tests['Acrylates/Methacrylates (Acute toxicity)']=humanBuiltQuery(create_test(), acry_words)
def create_test():
    acrylate=Chem.MolFromSmarts('[CH2]=[CH]C(=O)O')
    methacrylate=Chem.MolFromSmarts('[CH2]=C([CH3])C(=O)O')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        logp=x['logp']
        return (mol.HasSubstructMatch(acrylate) or mol.HasSubstructMatch(methacrylate)) and logp>5 and logp<8 and mw<1000
    return test
acry_words2 = {'qtype':'LogicalQuery', 'logic':'And', \
              'subqueries': [{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},\
                             {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'GreaterThan', 'value':5}, \
                             {'qtype':'b:ParameterQuery', 'prop':'Log Kow', 'operand':'LessThan', 'value':8}, \
                             {'qtype':'LogicalQuery', 'logic':'Or',\
                              'subqueries':[{'qtype':'b:StructureQuery', 'smart':'[CH2]=[CH]C(=O)O'},\
                                            {'qtype':'b:StructureQuery', 'smart':'[CH2]=C([CH3])C(=O)O'}]}]}
all_tests['Acrylates/Methacrylates (Chronic toxicity)']=humanBuiltQuery(create_test(), acry_words2)

#Epoxides
#Grace's advice
def create_test():
    strained_ring =Chem.MolFromSmarts('C1[O,N]C1')
    def test(x):
        mol=x['mol']
        mw=x['mol_weight']
        return mw<1000 and (mol.HasSubstructMatch(strained_ring) )
    return test
epox_words = {'qtype':'LogicalQuery', 'logic':'And',\
              'subqueries':[{'qtype':'b:ParameterQuery', 'prop':'Molecular Weight', 'operand':'LessThan', 'value':1000},\
                            {'qtype':'LogicalQuery', 'logic':'Or',\
                             'subqueries':[{'qtype':'b:StructureQuery', 'smart':'C1OC1'},\
                                           {'qtype':'b:StructureQuery', 'smart':'C1CN1'}]}]}
all_tests['Epoxides']=humanBuiltQuery(create_test(), epox_words)


# It looks like George never finished building this category, so for now it is being deprecated. 
del all_tests['Persistent, Bioaccumulative and Toxic (PBT) Chemicals']

# Type handling for chemical inputs
def normalizeChemicals(chemicals):
    """ This function checks the input type for the given chemicals and outputs the same information as 
    a DataFrame for smooth handling during querying."""
    # Check and handle input data types
    if type(chemicals) is dict:
        chemicals = pd.DataFrame(chemicals, index = [chemicals['dsstox_sid']])
    elif type(chemicals) is list:
        build_a_df = pd.DataFrame(chemicals[0], index = [chemicals[0]['dsstox_sid']])
        for chem in chemicals[1:]:
            chem = pd.DataFrame(chem, index = [chem['dsstox_sid']])
            build_a_df = pd.concat([build_a_df, chem])
        chemicals = build_a_df
    elif isinstance(chemicals, pd.DataFrame):
        pass
    else:
        raise TypeError("Chemicals must be supplied as a DataFrame, Dictionary, or list of Dictionaries")
    return chemicals

#Error Handling for all queries
def checkForAttributes(x):
    """This function ensures that all necessary attributes have been provided in the appropriate data types. Primarily
    checks that string values are not fully missing and that numerical columns contain numbers."""
    # Check inputs for all necessary attributes and data
    attributes = {'logp':'float64', 'ws':'float64', 'mol_weight':'float64', 'smiles':'object', 'dsstox_sid':'object', 'mol':'object'}
    key_counter = 0
    type_counter = 0
    bad_list = []
    for attr in attributes:
        try:
            x[attr]
        except KeyError:
            print(f"KeyError: {attr} must be provided")
            key_counter += 1
            bad_list.append(attr)
    for attr in bad_list:
        del attributes[attr]
    for attr in attributes:
        if x[attr].dtype != attributes[attr]:
            print(f"TypeError: {attr} must be provided as {attributes[attr]}. Please adjust the input accordingly.")
            type_counter += 1
    if key_counter > 0:
        raise KeyError('One or more attributes is missing. See printed statement(s).')
    elif type_counter > 0:
        raise TypeError('One or more attributes is of the wrong type. See printed statement(s).')
    else: 
        pass
    return None

def queryAll(chemicals, boolean_outputs = False):
    """This function will query every category for every chemical of interest.
       Inputs:
    - chemicals: A DataFrame, Dictionary, or list of Dictionaries of Chemicals and their attributes, including dsstox_sid, 
     smiles, logp, ws, mol_weight, and RDKIT MolfromSmiles (labelled as 'mol'). There must be keys or column names
     to match each these attribute titles.
      
       Returns:
    - category_df: A DataFrame of chemicals and their category memberships"""
    # Manage input type of chemicals
    chemicals = normalizeChemicals(chemicals)
    
    #Manage column types in chemicals DataFrame and ensure all attributes are present
    checkForAttributes(chemicals)

    categories = all_tests.keys()
    chem_dict = {'chemicals':chemicals['dsstox_sid']}
    for category in categories:
        chem_dict[category] = []
        for _, info in chemicals.iterrows():
            try:
                chem_dict[category].append(all_tests[category].query(info))
            except TypeError:
                print(f"{category} are still a problem")
                chem_dict[category].append(-1)
    category_df = pd.DataFrame(chem_dict, columns = chem_dict.keys())
    category_df = category_df.reset_index(drop = True)

    if boolean_outputs == True:
        pass
    else:
        for column in category_df.columns[1:]:
            category_df[column] = [int(i) for i in category_df[column]]

    return category_df

def listCategories(one_chem):
    """ Given a single chemical, this function will return a list of all categories to which that 
     chemical belongs. 
      
       Input: 
        - one_chem: A DataFrame or Dictionary representing a single chemical and its attributes, including dsstox_sid, 
     smiles, logp, ws, mol_weight, and RDKIT MolfromSmiles (labelled as 'mol'). There must be keys or column names
     to match each these attribute titles. 
      
       Output:
        - all_cats: A list of all categories to which the chemical belongs. """
    one_chem = normalizeChemicals(one_chem)
    checkForAttributes(one_chem)
    all_cats = []
    categories = all_tests.keys()
    for category in categories:
        if all_tests[category].query(dict(one_chem.iloc[0])):
            all_cats.append(category)
    return all_cats

possible_categories = list(all_tests.keys())

#These two functions simplify the necessary inputs for single queries and print trees

def singleQuery(category_title, one_chem):
    """This function takes in a category title and a single chemical entry and outputs whether the chemical is a 
    member of the specified category. 
    
           Inputs: 
        - category_title: A string representing a category title covered by these tests. Full list of acceptable
    categories is provided in the README.
        - one_chem: A DataFrame or Dictionary representing a single chemical and its attributes, including dsstox_sid, 
     smiles, logp, ws, mol_weight, and RDKIT MolfromSmiles (labelled as 'mol'). There must be keys or column names
     to match each these attribute titles. 
      
       Output:
        - boolean: A True/False value indicating the truth of the statement "one_chem is classified as category_title."

    At its core, singleQuery handles non-dictionary input types such as DataFrames and adds error handling for the print_tree feature of a Query
    instance and then runs the function."""
    one_chem = normalizeChemicals(one_chem)
    checkForAttributes(one_chem)
    if category_title not in possible_categories:
        print(f"{category_title} is not a valid category")
    return(all_tests[category_title].query(dict(one_chem.iloc[0])))

def printTree(category_title, one_chem = None, printer = True):
    """Given a category title and chemical, output the results of applying each query on the chemical
        to the console. For queries with subqueries, the subqueries will be disaplayed below the query in 
        indented lists. Can also be used to view query conditions without a chemical input. This is either a printed
        result OR a stored variable, depending on whether printer is set to True or False.
        
       Inputs: 
        - category_title: A string representing a category title covered by these tests. Full list of acceptable
        categories is provided in the README.
        - one_chem: Default None. A DataFrame or Dictionary representing a single chemical and its attributes, including dsstox_sid, 
        smiles, logp, ws, mol_weight, and RDKIT MolfromSmiles (labelled as 'mol'). There must be keys or column names
        to match each these attribute titles. 
        - printer: Default True. Whether or not you actually want to print the logic tree. With printer = False, this
        can output the tree as a savable string rather than actually printing it. 
      
       Output:
        - logic tree: A printed logic tree showing how classification decisions for the given category are made. If
        one_chem is provided, each branch will show whether one_chem satisfied the requirement or not for queries
        originating in the XML*. 

        *Queries that required hard-coded repair can be printed but will not show branch-by-branch results for application
        of constraints to one_chem.

        """
    if category_title not in possible_categories:
        print(f"{category_title} is not a valid category")
    if type(one_chem) != type(None):
        one_chem = normalizeChemicals(one_chem)
        checkForAttributes(one_chem)
        return(all_tests[category_title].print_tree(dict(one_chem.iloc[0]), printer = printer))
    else: 
        return(all_tests[category_title].print_tree(None, printer = printer))