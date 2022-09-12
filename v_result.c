
void v_result( const double *A, const double *B, double *C )
/**********************************************************/
/*    Vektorbetrag: C =  Vektor(B)-Vektor(A) == Vector(AB)*/
/**********************************************************/
{
         C[0]=B[0]-A[0];
         C[1]=B[1]-A[1];
         C[2]=B[2]-A[2];
}

