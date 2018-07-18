(* robotdistribution from parameters, reset distribution, and observations *)

ho[v_,m_,pars_]:=LogisticSigmoid[pars.{1,(2*(m-1)-200)/200,(2*(v-1)-200)/200}];

robotdistributiono[pf_,apars_,obs_]:= Block[
  {n=Length[obs],stub=Exp[apars[[1]]],pars=apars[[2;;]],cumf,a,a1,a2,b,b1,c},

  a=Table[Null,{n+1},{40}];
  
  cumf=Table[0,{n+1},{40}];
  Do[cumf[[i+1,obs[[i]]]]=1,{i,n}];
  cumf=Accumulate[cumf];

  a1=1;
  c={1};
  b={pf};

  a[[1]] = (c.b)/a1;


  Do[
    res=obs[[m-1]];
    a2=a1;
    a1=a[[m-1,res]];
    b1=b[[;;,res]];

    ht=Table[ho[v,m,pars],{v,m-1}];
    c=Prepend[(1-ht)*b1*c, (ht*b1).c]/a2;
    
    b= Table[(stub*pf + cumf[[m]]-cumf[[m-u]])/(stub+u) ,{u,0,m-1}];

    a[[m]] = (c.b)/a1;
    
   ,{m,2,n+1}];
a];
