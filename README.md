```
mkdir HEPData
cd HEPData/
curl -O https://repo.anaconda.com/archive/Anaconda3-2020.11-MacOSX-x86_64.sh
bash Anaconda3-2020.11-MacOSX-x86_64.sh <<< $'\nyes\n$PWD/anaconda3/\nyes\n'
rm Anaconda3-2020.11-MacOSX-x86_64.sh
source $PWD/anaconda3/bin/activate
rm ~/.bash_profile

conda update -n base -c defaults conda <<< $'y\n'
conda create -n py python=3 anaconda <<< $'y\n'
conda activate py
pip install hepdata_lib

hadd higgsCombineCombo.AsymptoticLimits.merged.MODELRPV.root `ls higgsCombineCombo.Asympto* | sort --version-sort -f`

python makeLimitHEPData.py
```
