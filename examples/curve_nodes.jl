Nnodes = 4
H = Vector{NodeFrame}(undef, Nnodes)
H[1] = NodeFrame(VEC3(),RV3(0.0,pi/6.0,0.0))
H[2] = NodeFrame(VEC3(1.0,0.0,0.0),RV3(0.0,-pi/6.0,0.0))
H[3] = NodeFrame(VEC3(1.0,1.0,0.0),RV3(pi/6.0,-pi/6.0,0.0))
H[4] = NodeFrame(VEC3(1.0,1.0,1.0),RV3(-pi/6.0,pi/3.0,pi/6.0))