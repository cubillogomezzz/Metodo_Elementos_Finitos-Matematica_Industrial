function [Cxl Cyl] = calcular_matrices_stokes(myPDE);

tol=1e-10;

Nn = length(myPDE.Mesh.Nodes);
Ne = length(myPDE.Mesh.Elements);
Nnl = max(max(myPDE.Mesh.Elements(1:3,:)));

%%Rg = zeros(6);
%%RR = sparse(Nn,Nn);
Cxg = zeros(6);
Cx = sparse(Nn,Nn);
Cyg = zeros(6);
Cy = sparse(Nn,Nn);
Cxg2 =zeros(6,3);
Cxl = sparse(Nn,Nnl);
Cyg2 =zeros(6,3);
Cyl = sparse(Nn,Nnl);

g=zeros(7,3);
g(1,:)=2;
g(6,:)=1;
g(7,:)=0;

for el = 1:Ne
	elemento=myPDE.Mesh.Elements(1:3,el);
	g(2,:)=myPDE.Mesh.Nodes(1,elemento)';
	g(4,:)=myPDE.Mesh.Nodes(2,elemento)';
	g(3,:)=g(2,[2 3 1]);
	g(5,:)=g(4,[2 3 1]);

	v1 = g([2 4],2) - g([2 4],1);
	v2 = g([2 4],3) - g([2 4],1);
	area = det([v1 v2])*0.5;

	hmin = 2*max(norm(v1),norm(v2));

	un_elem = createpde;
	geometryFromEdges(un_elem,g);
	generateMesh(un_elem,'hmin',hmin,'geometricorder','quadratic');

	aa = [];
	for nn=1:6
		phi=zeros(6,1);phi(nn)=1;results = createPDEResults(un_elem,phi);

		aa = [aa results];
	end

	xx = un_elem.Mesh.Nodes(1,1:6);
	yy = un_elem.Mesh.Nodes(2,1:6);
	area = det([xx(2)-xx(1) xx(3)-xx(1); yy(2)-yy(1) yy(3)-yy(1)])*0.5;

	for i=1:6
		for j=1:6
%%			Rg(j,i) = (aa(i).XGradients(4:6)'*aa(j).XGradients(4:6) + aa(i).YGradients(4:6)'*aa(j).YGradients(4:6))*area/3;
			Cxg(j,i) = aa(i).NodalSolution(4:6)'*aa(j).XGradients(4:6)*area/3;
			Cyg(j,i) = aa(i).NodalSolution(4:6)'*aa(j).YGradients(4:6)*area/3;
		end
	end

	aux = myPDE.Mesh.Elements(:,el);
	aux2 = un_elem.Mesh.Elements(:,1);


	ind1 = find(abs(myPDE.Mesh.Nodes(1,aux(1)) - un_elem.Mesh.Nodes(1,:)) < tol & abs(myPDE.Mesh.Nodes(2,aux(1)) - un_elem.Mesh.Nodes(2,:)) < tol);
	ind2 = find(abs(myPDE.Mesh.Nodes(1,aux(2)) - un_elem.Mesh.Nodes(1,:)) < tol & abs(myPDE.Mesh.Nodes(2,aux(2)) - un_elem.Mesh.Nodes(2,:)) < tol);
	ind3 = find(abs(myPDE.Mesh.Nodes(1,aux(3)) - un_elem.Mesh.Nodes(1,:)) < tol & abs(myPDE.Mesh.Nodes(2,aux(3)) - un_elem.Mesh.Nodes(2,:)) < tol);
	ind4 = find(abs(myPDE.Mesh.Nodes(1,aux(4)) - un_elem.Mesh.Nodes(1,:)) < tol & abs(myPDE.Mesh.Nodes(2,aux(4)) - un_elem.Mesh.Nodes(2,:)) < tol);
	ind5 = find(abs(myPDE.Mesh.Nodes(1,aux(5)) - un_elem.Mesh.Nodes(1,:)) < tol & abs(myPDE.Mesh.Nodes(2,aux(5)) - un_elem.Mesh.Nodes(2,:)) < tol);
	ind6 = find(abs(myPDE.Mesh.Nodes(1,aux(6)) - un_elem.Mesh.Nodes(1,:)) < tol & abs(myPDE.Mesh.Nodes(2,aux(6)) - un_elem.Mesh.Nodes(2,:)) < tol);


	ind = [ind1 ind2 ind3 ind4 ind5 ind6];
	Cxg = Cxg(ind,ind);
	Cyg = Cyg(ind,ind);

	Cxg2(:,1) = Cxg(:,1) + .5*Cxg(:,4) + .5*Cxg(:,6);
	Cxg2(:,2) = Cxg(:,2) + .5*Cxg(:,4) + .5*Cxg(:,5);
	Cxg2(:,3) = Cxg(:,3) + .5*Cxg(:,5) + .5*Cxg(:,6);

	Cyg2(:,1) = Cyg(:,1) + .5*Cyg(:,4) + .5*Cyg(:,6);
	Cyg2(:,2) = Cyg(:,2) + .5*Cyg(:,4) + .5*Cyg(:,5);
	Cyg2(:,3) = Cyg(:,3) + .5*Cyg(:,5) + .5*Cyg(:,6);

	Cxl(aux,aux(1:3)) = Cxl(aux,aux(1:3)) + Cxg2;
	Cyl(aux,aux(1:3)) = Cyl(aux,aux(1:3)) + Cyg2;
end

