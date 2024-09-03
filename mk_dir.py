from pathlib import Path
from decimal import Decimal, getcontext
import matplotlib.pyplot as plt

import numpy as np
import pandas as pd
import sys

#This script creates directories and conf files for mc



if (len(sys.argv)!=2):
    print("wrong number of arguments")
    exit(2)
rowNum=0#int(sys.argv[1])
unitCellNum=int(sys.argv[1])

inParamFileName="./V_1d_inv_12_6_2minParams.csv"



inDf=pd.read_csv(inParamFileName)
oneRow=inDf.iloc[rowNum,:]
a=float(oneRow.loc["a"])
b=float(oneRow.loc["b"])
A1=float(oneRow.loc["A1"])
A2=float(oneRow.loc["A2"])

diff1=float(oneRow.loc["diff1"])
diff2=float(oneRow.loc["diff2"])
sigma1=float(oneRow.loc["sigma1"])
sigma2=float(oneRow.loc["sigma2"])
f=float(oneRow.loc["f"])
d2=float(oneRow.loc["d2"])

d1=(2*a/b)**(1/6)

a1=d1-diff1
a2=d1+diff2

def V1(x):
    val1=a/x**12-b/x**6

    val2=-A1*np.exp(-1/sigma1**2*(x-a1)**2)-A2*np.exp(-1/sigma2**2*(x-a2)**2)


    val=val1+val2

    return val
# x=1
#
# print(V1(x))
# print("d1="+str(d1))
# print("a1="+str(a1))
# print("a2="+str(a2))
# xValsAll=np.linspace(a1,a2*1.5,100)
# yValsAll=[V1(x) for x in xValsAll]
#
# plt.figure()
# plt.plot(xValsAll,yValsAll,color="black")
# plt.savefig("tmp.png")
# plt.close()
# print(d2)

def format_using_decimal(value, precision=10):
    # Set the precision higher to ensure correct conversion
    getcontext().prec = precision + 2
    # Convert the float to a Decimal with exact precision
    decimal_value = Decimal(str(value))
    # Normalize to remove trailing zeros
    formatted_value = decimal_value.quantize(Decimal(1)) if decimal_value == decimal_value.to_integral() else decimal_value.normalize()
    return str(formatted_value)
TVals=[0.1,0.3,0.5,1,2,3,4,5,6,7,8,9,12,17,22,27,32]
# TVals=[3]
print("unitCellNum="+str(unitCellNum))
dataRoot="./dataAll/"
dataOutDir=dataRoot+"/dataAllUnitCell"+str(unitCellNum)+"/row"+str(rowNum)+"/"

TDirsAll=[]
TStrAll=[]
# print(TDirsAll)
for k in range(0,len(TVals)):
    T=TVals[k]
    # print(T)

    TStr=str(T)#format_using_decimal(T)
    TStrAll.append(TStr)
    TDir=dataOutDir+"/T"+TStr+"/"
    TDirsAll.append(TDir)
    Path(TDir).mkdir(exist_ok=True,parents=True)



def contents_to_conf(k):
    """

    :param k: index of T
    :return:
    """

    contents=[
        "#This is the configuration file for mc computations\n",
        "\n"
        "potential_function_name=V_1d_inv_12_6_2min\n",
        "\n" ,
        "#parameters of coefficients\n",
        "#parameter_row=row0\n",
        "#the following row is the values of a, b, A1, A2, a1, a2, sigma1, sigma2,f, d2\n"
        "coefs=["+format_using_decimal(a)+","+format_using_decimal(b)+","\
        +format_using_decimal(A1)+", "+format_using_decimal(A2)+", "+format_using_decimal(a1)\
        +", "+format_using_decimal(a2)+", "+format_using_decimal(sigma1)+", "+str(sigma2)\
        +", "+format_using_decimal(f)+", "+format(d2)+"]\n",
        "\n",
        "#Temperature\n",
        "T="+TStrAll[k]+"\n",
        "\n",
        "#unit cell number\n",
        "unitCellNum="+str(unitCellNum)+"\n",
        "\n",
        "erase_data_if_exist=False\n",
        "\n",
        "search_and_read_summary_file=True\n"
        "\n",
        "#For the observable name, only digits 0-9, letters a-zA-Z, underscore _ are allowed\n",
        "\n",
        "observable_name=U_dist\n",
        "\n",
        "effective_data_num_required=1000\n",
        "\n",
        "sweep_to_write=1000000\n",
        "\n",
        "#within each flush,  sweep_to_write mc computations are executed\n",
        "\n",
        "default_flush_num=15\n",
        "\n",
        "h=5e-2\n",
        "\n",
        "sweep_multiple=10\n"



    ]

    outConfName=TDirsAll[k]+"/run_T"+TStrAll[k]+".mc.conf"
    with open(outConfName,"w+") as fptr:
        fptr.writelines(contents)



for k in range(0,len(TDirsAll)):
    contents_to_conf(k)