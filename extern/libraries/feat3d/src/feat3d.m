      FEATFLOW=/home/cerberus5/featflow/sources/featflow1.2_new/featflow

         LIB=$FEATFLOW/object/libraries/libultra2/libfeat3d.a
      COMOPT="-xO5 -xtarget=ultra2 -dalign -xlibmil -fsimple=2 -Bstatic -depend -xlibmopt -lmvec -lcx -xarch=v8plusa -xsafe=mem -xcache=16/32/1:1024/64/1" 

for i in  *.f
do
      f77 -c $COMOPT -o o/$i.o $i
      ar rv $LIB o/$i.o 
      rm o/$i.o
      echo
done
      ranlib $LIB
