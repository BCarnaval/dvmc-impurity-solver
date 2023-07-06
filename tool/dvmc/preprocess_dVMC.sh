#!/usr/bin/env bash

# mVMC preparations :
dvmcdry.out StdFace.def

date
sed -i "s/NExUpdatePath  0/NExUpdatePath  1/"     modpara.def
sed -i "s/NDataQtySmp    1/NDataQtySmp    20/"    modpara.def

cp namelist.def namelist_G.def
cp modpara.def  modpara_G.def
cp namelist.def namelist_CA.def
cp modpara.def  modpara_CA.def
sed -i "s/CDataFileHead  zvo/CDataFileHead  zbo/"           modpara.def
sed -i "s/NVMCCalMode    0/NVMCCalMode    3/"               modpara_G.def
sed -i "s/ModPara  modpara.def/ModPara  modpara_G.def/"     namelist_G.def
sed -i "s/NVMCCalMode    0/NVMCCalMode    1/"               modpara_CA.def
sed -i "s/ModPara  modpara.def/ModPara  modpara_CA.def/"    namelist_CA.def


echo "      Excitation  excitation.def" >> namelist_G.def
makeExcitation.py
