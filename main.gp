encodegln(s,n)={
  my(v);
  v=[if(x==32,0,x-96)|x<-Vec(Vecsmall(s))];
  if(#v>n^2,warning("string truncated to length ",n^2));
  v = Vec(v,n^2);
  matrix(n,n,i,j,v[(i-1)*n+j]);
}

decodegln(m, n)={
  V = vector(n * n, i, 0);
  k = 1;
  for(i = 1, n, for(j = 1, n, V[k] = if(m[i,j] == 0, 32, m[i,j] + 96); k += 1));  
  Strchr(vecextract(V, "1..143"));
}

\\ algo d'exponentiation rapide itératif
fast_exp(m, n) =  {
	if (n == 0, m);
	acc = matid(matsize(m)[1]);
	while (n > 1, if (n%2 == 0, m = m^2; n = n / 2, acc = m * acc; m = m^2; n = (n - 1) / 2));
	acc = m * acc;
	acc;
};

\\ Idée: on exploite le fait que G=GL_12(Z/27Z) étant un groupe fini (non abélien),
\\ si k entier est premier avec l'ordre de ce groupe |G| alors l'application g|->g^k
\\ est une bijection sur G (découle directement de Bezout, g|->g^k définit une injection
\\ donc une bijection car G est fini).
\\ En particulier, il existe une unique racine k-ième dans G pour chaque k premier avec |G|.
\\ Calcul de la racine 65537 d'un matrice M (pour garantir l'unicité il faut que l'ordre
\\ du groupe soit premier avec 65537, ce qui est le cas), par le biais de la méthode de Newton :
mat_root(M,k)={ \\ k: nombre d'itérations
	my(X,i);
	X=matid(12);
	i=0;
	while(i<k, X=(X - (1/65537) * (fast_exp(X,65537) - M) * fast_exp(X^(-1), 65536)); i+=1);
	return(X);
}
\\ Il s'avère que si l'on se place dans Z/27Z alors la matrice n'est pas forcément inversible.
\\ Si jamais des coefficients dans F27 étaient employés pour l'exponentiation de la matrice d'origine,
\\ cette méthode pourrait fonctionner.
\\ Cette méthode ne fonctionne pas.


\\ Détermine le plus petit entier n pour lequel M^n=Id
ordreIdempotence(M, p) = {
    \\ L'ordre du groupe est product(27^12 - 27^k , k = 0 .. 11) = 125710791285604678382440447016150453040087085604527013876459165108951082279749262734927264400156938203820131742120111369101530910473068953664344868516622552271864320618748133936132438224326079460048633856000
    \\ l'élément M du groupe divise l'ordre du groupe,
    \\ on pourrait élever la matrice à la puissance chaque diviseur de l'ordre de ce groupe, bien que très lent
    R=Mod(M,p);
    k = 1;
    id = Mod(matid(12), p);
    while(id != R, R *= M; k+=1);
    return(k);
}


ciphertext=readstr("input.txt")[1];

C=Mod(encodegln(ciphertext,12), 27);



\\ On obtient C^o=(M^65537)^o=M^(65537*o)=Id
\\ De la relation de Bezout 65537u+ov=1 => 65537u=1-ov
\\ On en tire : C^u = M^(65537*u)=M^(1-ov)=M.M^(ov)=M.Id=M
o = ordreIdempotence(C, 27);
u=bezout(65537, o)[1];
C = fast_exp(C, u);
print(decodegln(lift(C),12));
