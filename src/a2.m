AA0=1/Sqrt[1+x^2+y^2]
AA1=x/Sqrt[1+x^2+y^2]
AA2=y/Sqrt[1+x^2+y^2]

AA={{AA2},{AA1},{AA0}}

m={{Sqrt[6/10], -Sqrt[2/6], Sqrt[1/15]}, {Sqrt[3/10], Sqrt[1/6], -Sqrt[8/15]}, {Sqrt[1/10], Sqrt[3/6], Sqrt[6/15]}}
im=Inverse[m]
a=im.AA
a/.{x->Sqrt[3],y->Sqrt[6]}
a/.{x->2.04,y->2.93}
a/.{x->2.44,y->3.53}
a/.{x->1.64,y->2.33}

Print[im.AA]

aaa={{a1},{a2},{a3}}
A=m.aaa

Print[A]


A/.{a1->1,a2->0,a3->0}
A/.{a1->0,a2->1,a3->0}
A/.{a1->0,a2->0,a3->1}
