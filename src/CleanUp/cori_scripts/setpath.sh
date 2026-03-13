module load cray-petsc
module load gsl
MODELGPATH=/global/homes/d/dteaney/langevinA/ModelG
MODELGPATHPY=/global/homes/d/dteaney/langevinA/python
O4MODELPATH=/global/homes/d/dteaney/o4/src
echo $MODELGPATH
echo $O4MODELPATH
export MODELGPATH
export O4MODELPATH
export PYTHONPATH=$MODELGPATH:$O4MODELPATH:$O4MODELGPATHPY


