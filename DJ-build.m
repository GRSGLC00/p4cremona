// ==========================================
// STEP 1: Define Rings and Setup
// ==========================================
Q := Rationals();
R2<x,y,z> := PolynomialRing(Q, 3);
R4<y0,y1,y2,y3,y4> := PolynomialRing(Q, 5);

// Define the points
e1:=[0,1,0]; e2:=[0,0,1]; e3:=[1,1,1]; e4:=[1,2,3]; e5:=[1,-1,2];
q:=[2,1,2];

// The generic conic C: d*x^2 + e*x*y + f*x*z + g*y*z = 0 (no y^2, z^2 since e1, e2 are on it)
// We find the coefficients passing through e3, e4, e5;
// Solving gives the conic:
C := 7*x^2 - 11*x*y + x*z + 3*y*z;

// The 5 cubics defining the anticanonical map
Y0 := x*C;
Y1 := y*C;
Y2 := z*C;
// We need two more general cubics passing through e1..e5
Y3 := 4*x^2*y + x*y^2 - 5*x^2*z; 
Y4 := 3*x^2*z + y*z^2 - 4*x*y*z;
Y_seq := [Y0, Y1, Y2, Y3, Y4];

// ==========================================
// STEP 2: The 12 Monomials in P4 containing L
// ==========================================
all_quads := MonomialsOfDegree(R4, 2);
bad_quads := [y3^2, y3*y4, y4^2];
M := [m : m in all_quads | not m in bad_quads]; 

// ==========================================
// STEP 3: Pull back to P2 and divide by C
// ==========================================
phi := hom< R4 -> R2 | Y_seq >;
V := [ ExactQuotient(phi(m), C) : m in M ];

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
NS := Nullspace(Transpose(Mat));

print "=== 5-DIMENSIONAL SYSTEM W ===";
print "Dimension of the System W (should be 5):", Dimension(NS);


W_quadrics := [];
i:=0;
for vec in Basis(NS) do
    quadric := &+[ vec[i] * M[i] : i in [1..12] ];
    Append(~W_quadrics, quadric);
    printf "W_%o = %o\n", i, quadric;
    i:= i+1;
end for;