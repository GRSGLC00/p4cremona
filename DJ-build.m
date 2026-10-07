// ==========================================
// STEP 1: Define Rings and Setup
// ==========================================
Q := Rationals();
R2<x,y,z> := PolynomialRing(Q, 3);
R4<[u]> := ProjectiveSpace(Q, 4);
R := CoordinateRing(R4);

// Define the points
pts := [
    [0,1,0],
    [0,0,1],
    [1,1,1],
    [1,2,3],
    [1,-2,5]
];
q:=[2,1,3];
monomials := [x^2, y^2, z^2, x*y, x*z, y*z];
eval_matrix := Matrix(Q, [ [Evaluate(m, pt) : m in monomials] : pt in pts ]);
N := Nullspace(Transpose(eval_matrix));
if Dimension(N) eq 0 then
    error "No conic passes through all given points.";
end if;
// Reconstruct the polynomial in R from the nullspace basis vector
v := Basis(N)[1];
C := &+[ v[i] * monomials[i] : i in [1..#monomials] ];
// The 5 cubics defining the anticanonical map
monomials_deg3 := [x^3, x^2*y, x^2*z, x*y^2, x*y*z, x*z^2, y^3, y^2*z, y*z^2, z^3];
eval_matrix3 := Matrix(Q, [ [Evaluate(m, pt) : m in monomials_deg3] : pt in pts ]);
// W is the 5-dimensional space of cubics vanishing at the 5 points
W := Nullspace(Transpose(eval_matrix3));
Y_seq := [];
Y_seq[1]:=x*C;
Y_seq[2]:=y*C;
Y_seq[3]:=z*C;
cubic1:=Basis(W)[1];
cubic2:=Basis(W)[2];
Y_seq[4]:=&+[ cubic1[i] * monomials_deg3[i] : i in [1..10] ];
Y_seq[5]:=&+[ cubic2[i] * monomials_deg3[i] : i in [1..10] ];

// ==========================================
// STEP 2: The 12 Monomials in P4 containing L
// ==========================================
all_quads := MonomialsOfDegree(R, 2);
bad_quads := [u[4]^2, u[4]*u[5], u[5]^2];
M := [m : m in all_quads | not m in bad_quads]; 

// ==========================================
// STEP 3: Pull back to P2 and divide by C
// ==========================================
phi := hom< R -> R2 | Y_seq >;
V := [ ExactQuotient(phi(m), C) : m in M ];
V;

// ==========================================
// STEP 4: Set up the Linear System 
// ==========================================
rows := []; 

// Condition 1: Pass through q
Append(~rows, [ Evaluate(v, q) : v in V ]);

// Conditions 2-7: Triple point at p = [1,0,0]
// We simply extract the coefficients of the 6 monomials that MUST be zero
forbidden_monomials := [x^4, x^3*y, x^3*z, x^2*y^2, x^2*y*z, x^2*z^2];

for mon in forbidden_monomials do
    row := [ MonomialCoefficient(v, mon) : v in V ];
    Append(~rows, row);
end for;

Mat := Matrix(Q, 7, 12, rows);

// ==========================================
// STEP 5: Reconstruct the Quadrics W
// ==========================================
// We compute the solution of the linear system, finding the five quadrics of P^4
NS := Nullspace(Transpose(Mat));

print "=== 5-DIMENSIONAL SYSTEM W ===";
print "Dimension of the System W:", Dimension(NS);


W_quadrics := [];
i:=0;
for vec in Basis(NS) do
    quadric := &+[ vec[i] * M[i] : i in [1..12] ];
    Append(~W_quadrics, quadric);
    printf "W_%o := %o;\n", i, quadric;
    i:= i+1;
end for;
