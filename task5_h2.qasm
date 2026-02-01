OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];

// Task 5: exp(i pi/4 (XX + YY + ZZ)) = exp(i pi/4) * SWAP
cx q[0], q[1];
cx q[1], q[0];
cx q[0], q[1];
