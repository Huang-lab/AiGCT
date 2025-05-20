import context  # noqa: F401

from aigct.model import VariantId  # noqa: F401
from dataclasses import dataclass, field
import os
import pandas as pd
from sklearn.metrics import (
    roc_auc_score, roc_curve, precision_recall_curve, auc)


print(type(VariantId))


df = pd.read_csv("temp/ASD_630.csv")



@dataclass
class TableDef:
    folder: str
    file_name: str
    pk_columns: list[str]
    non_pk_columns: list[str]
    columns: field(init=False) = None,
    full_file_name: field(init=False) = None

    def __post_init__(self):
        self.columns = self.pk_columns + self.non_pk_columns
        self.full_file_name = os.path.join(self.folder, self.file_name)


# TableDef("a", "b", ["a","b"], ["a","d"])      

class par:

    def __init__(self, prop):
        self._prop = prop

    @property
    def prop(self):
        return self._prop
    

class pars(par):

    def __init__(self, prop):
        super().__init__(prop)

inst = pars('a')
inst.prop

di = {"a": "b", "c": "d"}
dr = {v: k for k, v in di.items()}


df = pd.DataFrame({"a": [1,1,2,2,3,4], "b": [1,2,3,4,5,6]})
dfg = df.groupby("a")

import numpy as np

df = pd.DataFrame({'A': ["a","a",np.nan,np.nan,"a"],
                    'B': [1, 2, 3, 4,1],
                    'C': [4, 6, 5, 5,4],
                    'D': ["a","c","b",np.nan,"a"]})
g1 = df.groupby('A', group_keys=False)
g2 = df.groupby('A', group_keys=True)


dfqr = df.query('A != D')

xx = df.drop_duplicates()

def lf(x):
    ret = x/x.sum()
    return ret

def lf1(x):
    m1 = x['B'].max() + x['C'].max()
    # s = pd.Series({"A": x['A'], "B": x['B'].max(), "C": x['C'].max(), "D": m1})
    s = pd.Series({"B": x['B'].max(), "C": x['C'].max(), "D": m1})
    return s

gg = g1[['B','C']]

# g1[['B', 'C']].apply(lambda x: x / x.sum())
#g1[['B', 'C']].apply(lf)
gr = g2[['B', 'C']].apply(lf1)

roc_curve([1,1],[.3,.4])

try:
    roc_auc_score([1,1],[.3,.4])
except ValueError as ve:
    a = str(ve)
except Exception as e:
    b = str(e)

try:
    prs, recs, thresholds = precision_recall_curve([1,1],[.2,.3])
    auc(recs, prs)
except ValueError as ve:
    d = str(ve)
except Exception as e:
    f = str(e)



with open("junk.lis", "w") as out:
    out.write("abc\n\n")
    out.write("xbc")

def tfunc(*args):
    for a in args:
        print(a)
    return [1,2,3]

x,y,z = tfunc(5,6,7)
x

import matplotlib.colors as mcolors

a = mcolors.CSS4_COLORS

print('start')
for r in df:
    print(r)

for i,r in df.iterrows():
    print(str(i) + ": " + str(r))

pass

import os
y = ("a","b")
x = os.path.join(*y)
pass

dic = [{"a":1,"b":2},{"a":5,"b":6}]

df3= pd.DataFrame([{"a":1,"b":2},{"a":5,"b":6}],columns=["A","B"])
pass

from aigct.repository import (
    VARIANT_PK_COLUMNS
)
v1 = pd.read_csv("data/variant.csv")
v1s = v1.copy()
vd = pd.read_csv("temp/variant.csv")
vd2 = vd.query("not ALLELE_FREQUENCY.isna()")

v1 = v1.set_index(VARIANT_PK_COLUMNS)
vd2 = vd2.set_index(VARIANT_PK_COLUMNS)[
    ["ALLELE_FREQUENCY_SOURCE", "ALLELE_FREQUENCY"]]
v1.update(vd2)

v1 = v1.reset_index()

v1.to_csv("data/variant.csv",index=False)


df = pd.DataFrame({
    'group': ['A', 'A', 'B', 'B', 'B', 'C'],
    'col1': [1, 1, 2, 2, 3, 4],
    'col2': ['x', 'y', 'x', 'x', 'z', 'z'],
    'col3': [True, False, True, True, False, True]
})

# Method 1: Group by one column and get unique counts for specific columns
result = df.groupby('group').apply(
    lambda x: x[['col1', 'col2', 'col3']].drop_duplicates().shape[0]
)

pass