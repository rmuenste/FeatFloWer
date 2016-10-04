

set(src_feat3d
src/ab03.f
src/ab07.f
src/ab09.f
src/ab17.f
src/ap3.f
src/ap7.f
src/ap9.f
src/cb3h.f
src/e010.f
src/e011.f
src/e013.f
src/e030.f
src/e031.f
src/gorsm.f
src/gowsm.f
src/mp311.f
src/mp330.f
src/mp331.f
src/mr311.f
src/mr330.f
src/mr331.f
src/ndfg.f
src/ndfgl.f
src/ndfl.f 
src/orsc.f
src/sadb.f
src/sb0.f
src/sbael.f
src/sba.f
src/sbca.f
src/sbcb.f
src/sbeel.f
src/sbe.f 
src/sbvel.f
src/sbv.f
src/sedb.f
src/svc.f
src/vb0.f
src/xmab07.f
src/xmab17.f
src/xmap.f 
src/xmora.f
src/xmors.f
src/xmowa.f
src/xmows.f
src/xmsb0.f
src/xmscl.f
src/xors.f
src/xows.f
src/xsb0x.f
src/zfill.f
src/zvalue.f
)

# add the math library
if(FC_CUDA_SUPPORT)
  cuda_add_library(feat3d ${src_feat3d})
else(FC_CUDA_SUPPORT)
  add_library(feat3d ${src_feat3d})
endif(FC_CUDA_SUPPORT)
