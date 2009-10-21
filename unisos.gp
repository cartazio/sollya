/* Show timings on examples */

default(timer,1);

/* Need to increase the default a bit for some of the examples */

default(realprecision,100);

/* Return the "squared part" of a polynomial so that if
 *
 * s = squaredpart(p)
 *
 * then we have p = s^2 q where q is a polynomial without double roots.
 *
 */

squaredpart(p) =
{ local(n,rp,lp,fp);
  /* Suppose p = c * product of (x - a_k)^m_k */
  n = poldegree(p);
  /* Let rp[i] = c * product of (x - a_k)^{max(m_k-i,0)} */
  rp = listcreate(n);
  listput(rp,gcd(p,deriv(p)),1);
  for (i=2,n,listput(rp,gcd(rp[i-1],deriv(rp[i-1])),i));
  /* Let lp[i] = product of (x - a_k)^{1 if a_k >= i, 0 otherwise} */
  lp = listcreate(n+1);
  listput(lp,p/rp[1],1);
  for (i=2,n,listput(lp,rp[i-1]/rp[i],i));
  listput(lp,1,n+1);
  /* let fp[i] = product of (x - a_k)^{1 if a_k = i, 0 otherwise} */
  /* i.e. the product of factors with root multiplicity exactly i */
  fp = listcreate(n);
  for (i=1,n,listput(fp,lp[i]/lp[i+1],i));
  /* Now pick out maximal squarable factor from that */
  return(prod(i=1,n,fp[i]^(floor(i/2))))
}

/* Round at the 2^-i position */

myround(p,i) = round(2^i*p)/2^i;

/*
 * Main function
 */

sos(p) =
{
  local(s,c,q,n,m,t,e,r,k,ok,l,a,p_can,p_cnj,s1,s2,u,v,i,sqs,cfs);

  /* Factor out the repeated roots, separating squared part s    */
  /* Also get rid of the coefficient so the remaining q is monic */
  s = squaredpart(p);
  c = polcoeff(p,poldegree(p));
  q = p / (c * s^2);

  /* Let n = 2 m be the degree; fail if this is odd (can't be pos def) */
  /* If it's a constant (i.e. the polynomial is a constant multiple of */
  /* a perfect square) then return almost immediately.                 */
  n = poldegree(q);
  m = floor(n/2);
  if (2 * m != n,return);
  if (n == 0,
        print(p/s^2 " * (" s ")^2");
        return
  );

  /* let t = 1 + x^2 + ... + x^2m and r = q - e * t be safe perturbation */
  t = (x^(n+2)-1)/(x^2-1);
  e = 1;
  while((q-e*t != 0 & polsturm(q-e*t) != 0),e = e / 2); 
  e = e / 2;
  r = q - e * t;

  /* Choose rounding k of roots accurate enough to make the rest work */
  k = 0; ok = 0;
  until(ok,
    k = k + 1;
    l = pollead(r);
    a = polroots(r);

    p_can = prod(i = 1,poldegree(r)/2,x - a[2*i-1]);
    p_cnj = prod(i = 1,poldegree(r)/2,x - a[2*i]);
    s1 = myround(real((p_can + p_cnj) / 2),k);
    s2 = myround(real((p_can - p_cnj) / (2 * I)),k);
    u = r - l * (s1^2 + s2^2);
    v = e * t + u;

    ok = 1;
    for(i = 0,m,ok = ok *
       (polcoeff(v,2*i) >= abs(polcoeff(v,2*i+1))/4 + abs(polcoeff(v,2*i-1))))
  );

  /* Accumulate the final list of squares and coefficients */
  sqs = vector(2 * m + 3);
  cfs = vector(2 * m + 3);
  sqs[1] = s1; cfs[1] = l;
  sqs[2] = s2; cfs[2] = l;
  for (i = 0,m,
    sqs[i+3] = x^i;
    cfs[i+3] = polcoeff(v,2*i) -
               (abs(polcoeff(v,2*i+1))/4 + abs(polcoeff(v,2*i-1)))
  );
  for (i = 0,m-1,
    sqs[i+m+4] = x^i * (x + sign(polcoeff(v,2*i+1)) / 2);
    cfs[i+m+4] = abs(polcoeff(v,2*i+1))
  );

 /* Now put back in the factor from the initial decomposition */

  for (i = 1,2 * m + 3,
    sqs[i] = s * sqs[i];
    cfs[i] = c * cfs[i]
  );

 /* Sanity check */

  if (sum(i=1,2*m+3,cfs[i]*sqs[i]^2) - p != 0,print("Failure");return);
  for (i = 1,2*m+2,print(cfs[i] " * (" sqs[i] ")^2 +"));
  print(cfs[2*m+3] " * (" sqs[2*m+3] ")^2")
}

/*
 * Now generate a PSatz certificate over an interval directly
 */

psatz(p,a,b) =
{
  local(s,c,q,n,m,t,e,r,k,ok,l,ar,p_can,p_cnj,s1,s2,u,v,i,j,jj);
  local(d,d2,zesqs,ecfs,osqs,ocfs);

  /** Form range-reduced polynomial for the main computation **/
  d = poldegree(p);
  d2 = floor(d/2);
  pp = (1 + x^2)^d * subst(p,x,(a + b*x^2)/(1 + x^2));

  /* Factor out the repeated roots, separating squared part s    */
  /* Also get rid of the coefficient so the remaining q is monic */
  s = squaredpart(pp);
  c = polcoeff(pp,poldegree(pp));
  q = pp / (c * s^2);

  /* Let n = 2 m be the degree.                                        */
  n = poldegree(q);
  m = floor(n/2);

  /* Allocate vector for sum_i cfs[i] * sqs[i]^2 */
  sqs = vector(2 * m + 3);
  cfs = vector(2 * m + 3);

  /* If we just have a constant then it's rather trivial */
  if(n == 0,
      cfs[1] = q; sqs[1] = 1;
      cfs[2] = 0; sqs[2] = 0;
      cfs[3] = 0; sqs[3] = 0
  /* Otherwise use the main algorithm */
     ,
      /* let t = 1 + x^2 + ... + x^2m and r = q - e * t
         be a safe perturbation */
      t = (x^(n+2)-1)/(x^2-1);
      e = 1;
      while((q-e*t != 0 & polsturm(q-e*t) != 0),e = e / 2);
      e = e / 2;
      r = q - e * t;

      /* Choose rounding k of roots accurate enough to make the rest work */
      k = 0; ok = 0;
      until(ok,
        k = k + 1;
        l = pollead(r);
        ar = polroots(r);

        p_can = prod(i = 1,poldegree(r)/2,x - ar[2*i-1]);
        p_cnj = prod(i = 1,poldegree(r)/2,x - ar[2*i]);
        s1 = myround(real((p_can + p_cnj) / 2),k);
        s2 = myround(real((p_can - p_cnj) / (2 * I)),k);
        u = r - l * (s1^2 + s2^2);
        v = e * t + u;

        ok = 1;
        for(i = 0,m,ok = ok *
           (polcoeff(v,2*i) >=
                abs(polcoeff(v,2*i+1))/4 + abs(polcoeff(v,2*i-1))))
      );

      /* Accumulate the final list of squares and coefficients */
      sqs[1] = s1; cfs[1] = l;
      sqs[2] = s2; cfs[2] = l;
      for (i = 0,m,
        sqs[i+3] = x^i;
        cfs[i+3] = polcoeff(v,2*i) -
                   (abs(polcoeff(v,2*i+1))/4 + abs(polcoeff(v,2*i-1)))
      );
      for (i = 0,m-1,
        sqs[i+m+4] = x^i * (x + sign(polcoeff(v,2*i+1)) / 2);
        cfs[i+m+4] = abs(polcoeff(v,2*i+1))
      );
    );

  /* Now put back in the factor from the initial decomposition */
  for (i = 1,2 * m + 3,
    sqs[i] = s * sqs[i];
    cfs[i] = c * cfs[i]
  );

  /* Now split out the squares to give all even terms */
  esqs = vector(2 * m + 3);
  ecfs = vector(2 * m + 3);
  osqs = vector(2 * m + 3);
  ocfs = vector(2 * m + 3);

  for(i = 1,2 * m + 3,
        ecfs[i] = cfs[i];
        esqs[i] = 0;
        for(j = 0,poldegree(sqs[i]),
            esqs[i] = esqs[i] + x^j * polcoeff(sqs[i],2*j));
        ocfs[i] = cfs[i];
        osqs[i] = 0;
        for(j = 0,poldegree(sqs[i]),
            osqs[i] = osqs[i] + x^j * polcoeff(sqs[i],2*j+1));
  );

  /* Sanity check */
  if (sum(i=1,2*m+3,ecfs[i]*(subst(esqs[i],x,x^2)^2)) +
      x^2 * sum(i=1,2*m+3,ocfs[i]*(subst(osqs[i],x,x^2)^2)) -
      pp != 0,print("Failure: sanity check on intermediate poly");return);

  /* Now modify the terms to put them in terms of original variable */
  /* We need two cases according to whether the original degree is  */
  /* even or odd                                                    */

  for(i = 1,2 * m + 3,
        ecfs[i] = ecfs[i] / (b - a)^d;
        ocfs[i] = ocfs[i] / (b - a)^d;
        esqs[i] = (b - x)^d2 * subst(esqs[i],x,(x-a)/(b-x));
        if(2 * d2 == d,
           osqs[i] = (b - x)^(d2-1) * subst(osqs[i],x,(x-a)/(b-x)),
           osqs[i] = (b - x)^d2 * subst(osqs[i],x,(x-a)/(b-x))));

  /* Now another sanity check */

  if(2 * d2 == d,
     if(sum(i=1,2*m+3,ecfs[i]*esqs[i]^2) +
        (x - a) * (b - x) * sum(i=1,2*m+3,ocfs[i]*osqs[i]^2) != p,
        print("Failure: sanity check on final poly (even)"); return),
     if((b - x) * sum(i=1,2*m+3,ecfs[i]*esqs[i]^2) +
        (x - a) * sum(i=1,2*m+3,ocfs[i]*osqs[i]^2) != p,
        print("Failure: sanity check on final poly (odd)"); return));

  /* Now print out the final result */

  print("q = (");

  if(2 * d2 != d,print("(" b " - x) * ("),print("("));
  for (i = 1,2*m+2,print(ecfs[i] " * (" esqs[i] ")^2 +"));
  print(ecfs[2*m+3] " * (" esqs[2*m+3] ")^2) +");
  if(2 * d2 != d,print("(x - " a ") * ("),
                 print("(x - " a ") * (" b " - x) * ("));
  for (i = 1,2*m+2,print(ocfs[i] " * (" osqs[i] ")^2 +"));
  print(ocfs[2*m+3] " * (" osqs[2*m+3] ")^2)");

  print(");")

}



\r sosInput.gp

trap(,quit,psatz(p,a,b));

\q


