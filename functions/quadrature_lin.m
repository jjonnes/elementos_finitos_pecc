function Q = quadrature_lin(q)
quadrature_1_pts = [0.0  2.0];
quadrature_2_pts = [
    -0.577350269189626 1.0
     0.577350269189626 1.0];
quadrature_3_pts = [
    -0.774596669241483 5.0/9.0
     0.0               8.0/9.0
     0.774596669241483 5.0/9.0];
quadrature_4_pts = [
     -0.861136312 0.3478548
     -0.339981044 0.6521452
      0.339981044 0.6521452
      0.861136312 0.3478548];

if q == 1
    Q = quadrature_1_pts;
elseif q == 2
    Q = quadrature_2_pts;
elseif q == 3
    Q = quadrature_3_pts;
elseif q == 4
    Q = quadrature_4_pts;
else
    Q = 0; 
end
end