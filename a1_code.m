function [iter_arr,res_arr,A1,A2,A3,A4,A5,A6] = a1_code(Nx, Ny, Nz, sigma,f1,f2,f3,f4,f5,f6, thres)
%A_i: discretization matrix for F_i; f_i: boundary condition for F_i
dx = 2 / Nx; dy = 2 / Ny; dz = 2 / Nz;
x = linspace(-1+dx/2, 1-dx/2, Nx); %row vector
y = linspace(-1+dy/2, 1-dy/2, Ny);
z = linspace(-1+dz/2, 1-dz/2, Nz);
[X, Y, Z] = meshgrid(x,y,z);

A1 = zeros(Nx+2, Ny+2, Nz+2); A2 = A1; A3 = A1; A4 = A1; A5 = A1; A6 = A1;
n_iter = 0;

% Boundary conditions
[X_xy,Y_xy]=meshgrid(x,y); [X_xz,Z_xz]=meshgrid(x,z);[Y_yz,Z_yz]=meshgrid(y,z);
A1b = f1(Y_yz,Z_yz); A2b = f2(Y_yz,Z_yz); A3b = f3(X_xz,Z_xz); A4b = f4(X_xz,Z_xz); A5b=f5(X_xy,Y_xy); A6b=f6(X_xy,Y_xy);
A1(1,2:end-1,2:end-1) = A1b'; A2(end,2:end-1,2:end-1) = A2b';
A3(2:end-1,1,2:end-1) = A3b'; A4(2:end-1,end,2:end-1) = A4b';
A5(2:end-1,2:end-1,1) = A5b'; A6(2:end-1,2:end-1,end) = A6b';

% Residual
resA1 = -(A1(2:Nx+1,2:Ny+1,2:Nz+1) - A1(1:Nx,2:Ny+1,2:Nz+1)) / dx + sigma * ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A1(2:Nx+1,2:Ny+1,2:Nz+1));
resA2 = (A2(3:Nx+2,2:Ny+1,2:Nz+1) - A2(2:Nx+1,2:Ny+1,2:Nz+1)) / dx + sigma * ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A2(2:Nx+1,2:Ny+1,2:Nz+1));
resA3 = -(A3(2:Nx+1,2:Ny+1,2:Nz+1) - A3(2:Nx+1,1:Ny,2:Nz+1)) / dy + sigma * ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A3(2:Nx+1,2:Ny+1,2:Nz+1));
resA4 = (A4(2:Nx+1,3:Ny+2,2:Nz+1) - A4(2:Nx+1,2:Ny+1,2:Nz+1)) / dy + sigma * ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A4(2:Nx+1,2:Ny+1,2:Nz+1));
resA5 = -(A5(2:Nx+1,2:Ny+1,2:Nz+1) - A5(2:Nx+1,2:Ny+1,1:Nz)) / dz + sigma * ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A5(2:Nx+1,2:Ny+1,2:Nz+1));
resA6 = (A6(2:Nx+1,2:Ny+1,3:Nz+2) - A6(2:Nx+1,2:Ny+1,2:Nz+1)) / dz + sigma * ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A6(2:Nx+1,2:Ny+1,2:Nz+1));

res=norm(resA1, 'fro')^2 + norm(resA2, 'fro')^2 + norm(resA3, 'fro')^2 + norm(resA4, 'fro')^2 + norm(resA5, 'fro')^2 + norm(resA6, 'fro')^2;
res_arr = [];
iter_arr=[];

while res > thres && n_iter < 8000
    n_Newton_iter = 0;
    n_nonlin_eqs = 0;

    count=0;
    for i=2:Nx+1 
        for j=2:Ny+1
            for k=2:Nz+1
                count=count+1;
                resA1 = -(A1(i,j,k) - A1(i-1,j,k)) / dx + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A1(i,j,k)); 
                resA2 = (A2(i+1,j,k) - A2(i,j,k)) / dx + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A2(i,j,k));
                resA3 = -(A3(i,j,k) - A3(i,j-1,k)) / dy + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A3(i,j,k));
                resA4 = (A4(i,j+1,k) - A4(i,j,k)) / dy + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A4(i,j,k));
                resA5 = -(A5(i,j,k)-A5(i,j,k-1)) / dz + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A5(i,j,k));
                resA6 = (A6(i,j,k+1)-A6(i,j,k)) / dz + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A6(i,j,k));
                res = sqrt(resA1^2 + resA2^2 + resA3^2 + resA4^2+ resA5^2 + resA6^2);
                
                n_nonlin_eqs = n_nonlin_eqs + 1;
                
                while res > thres * 1e-3
                    Jac=[1/dx - sigma/6 + sigma, -sigma/6, -sigma/6, -sigma/6, -sigma/6, -sigma/6;
                        -sigma/6, 1/dx - sigma/6 + sigma, -sigma/6, -sigma/6, -sigma/6, -sigma/6;
                        -sigma/6, -sigma/6, 1/dy - sigma/6 + sigma, -sigma/6, -sigma/6, -sigma/6;
                        -sigma/6, -sigma/6, -sigma/6, 1/dy - sigma/6 + sigma, -sigma/6, -sigma/6;
                        -sigma/6, -sigma/6, -sigma/6, -sigma/6, 1/dz - sigma/6 + sigma, -sigma/6;
                        -sigma/6, -sigma/6, -sigma/6, -sigma/6, -sigma/6, 1/dz - sigma/6 + sigma];

                    dsol = Jac \ [resA1; resA2; resA3; resA4; resA5; resA6];
    
                    A1(i,j,k) = A1(i,j,k) + dsol(1);
                    A2(i,j,k) = A2(i,j,k) + dsol(2);
                    A3(i,j,k) = A3(i,j,k) + dsol(3);
                    A4(i,j,k) = A4(i,j,k) + dsol(4);
                    A5(i,j,k) = A5(i,j,k) + dsol(5);
                    A6(i,j,k) = A6(i,j,k) + dsol(6);
  
                 
                    resA1 = -(A1(i,j,k) - A1(i-1,j,k)) / dx + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A1(i,j,k)); 
                    resA2 = (A2(i+1,j,k) - A2(i,j,k)) / dx + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A2(i,j,k));
                    resA3 = -(A3(i,j,k) - A3(i,j-1,k)) / dy + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A3(i,j,k));
                    resA4 = (A4(i,j+1,k) - A4(i,j,k)) / dy + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A4(i,j,k));
                    resA5 = -(A5(i,j,k)-A5(i,j,k-1))/dz + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A5(i,j,k));
                    resA6 = (A6(i,j,k+1)-A6(i,j,k))/dz + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A6(i,j,k));
                    res = sqrt(resA1^2 + resA2^2 + resA3^2 + resA4^2+ resA5^2 + resA6^2);

                    

                    n_Newton_iter = n_Newton_iter + 1;
    
                    if (res > 1e+8)
                        error("The Newton iteration fails to converge.");
                    end
                end
            end
        end
    end

    for i=Nx+1:-1:2
        for j=Ny+1:-1:2
            for k=Nz+1:-1:2
                resA1 = -(A1(i,j,k) - A1(i-1,j,k)) / dx + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A1(i,j,k)); 
                resA2 = (A2(i+1,j,k) - A2(i,j,k)) / dx + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A2(i,j,k));
                resA3 = -(A3(i,j,k) - A3(i,j-1,k)) / dy + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A3(i,j,k));
                resA4 = (A4(i,j+1,k) - A4(i,j,k)) / dy + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A4(i,j,k));
                resA5 = -(A5(i,j,k)-A5(i,j,k-1))/dz + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A5(i,j,k));
                resA6 = (A6(i,j,k+1)-A6(i,j,k))/dz + sigma * ...
                    ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A6(i,j,k));
                res = sqrt(resA1^2 + resA2^2 + resA3^2 + resA4^2+ resA5^2 + resA6^2);


                n_nonlin_eqs = n_nonlin_eqs + 1;
                while res > 1e-3 * thres 
                    Jac=[1/dx - sigma/6 + sigma, -sigma/6, -sigma/6, -sigma/6, -sigma/6, -sigma/6;
                         -sigma/6, 1/dx - sigma/6 + sigma, -sigma/6, -sigma/6, -sigma/6, -sigma/6;
                         -sigma/6, -sigma/6, 1/dy - sigma/6 + sigma, -sigma/6, -sigma/6, -sigma/6;
                         -sigma/6, -sigma/6, -sigma/6, 1/dy - sigma/6 + sigma, -sigma/6, -sigma/6;
                         -sigma/6, -sigma/6, -sigma/6, -sigma/6, 1/dz - sigma/6 + sigma, -sigma/6;
                         -sigma/6, -sigma/6, -sigma/6, -sigma/6, -sigma/6, 1/dz - sigma/6 + sigma];
    
                    dsol = Jac \ [resA1; resA2; resA3; resA4; resA5; resA6];
        
                    A1(i,j,k) = A1(i,j,k) + dsol(1);
                    A2(i,j,k) = A2(i,j,k) + dsol(2);
                    A3(i,j,k) = A3(i,j,k) + dsol(3);
                    A4(i,j,k) = A4(i,j,k) + dsol(4);
                    A5(i,j,k) = A5(i,j,k) + dsol(5);
                    A6(i,j,k) = A6(i,j,k) + dsol(6);
    
         
                    resA1 = -(A1(i,j,k) - A1(i-1,j,k)) / dx + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A1(i,j,k)); 
                    resA2 = (A2(i+1,j,k) - A2(i,j,k)) / dx + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A2(i,j,k));
                    resA3 = -(A3(i,j,k) - A3(i,j-1,k)) / dy + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A3(i,j,k));
                    resA4 = (A4(i,j+1,k) - A4(i,j,k)) / dy + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A4(i,j,k));
                    resA5 = -(A5(i,j,k)-A5(i,j,k-1))/dz + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A5(i,j,k));
                    resA6 = (A6(i,j,k+1)-A6(i,j,k))/dz + sigma * ...
                        ((1/6)*(A1(i,j,k) + A2(i,j,k) + A3(i,j,k) + A4(i,j,k) + A5(i,j,k)+A6(i,j,k))-A6(i,j,k));
    
                    res = sqrt(resA1^2 + resA2^2 + resA3^2 + resA4^2+ resA5^2 + resA6^2);
    
                    n_Newton_iter = n_Newton_iter + 1;
    
                    if (res > 1e+8)
                        error("The Newton iteration fails to converge.");
                    end
                end
            end
        end
    end

    resA1 = -(A1(2:Nx+1,2:Ny+1,2:Nz+1) - A1(1:Nx,2:Ny+1,2:Nz+1)) / dx + sigma * ...
        ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A1(2:Nx+1,2:Ny+1,2:Nz+1));
    resA2 = (A2(3:Nx+2,2:Ny+1,2:Nz+1) - A2(2:Nx+1,2:Ny+1,2:Nz+1)) / dx + sigma * ...
        ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A2(2:Nx+1,2:Ny+1,2:Nz+1));
    resA3 = -(A3(2:Nx+1,2:Ny+1,2:Nz+1) - A3(2:Nx+1,1:Ny,2:Nz+1)) / dy + sigma * ...
        ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A3(2:Nx+1,2:Ny+1,2:Nz+1));
    resA4 = (A4(2:Nx+1,3:Ny+2,2:Nz+1) - A4(2:Nx+1,2:Ny+1,2:Nz+1)) / dy + sigma * ...
        ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A4(2:Nx+1,2:Ny+1,2:Nz+1));
    resA5 = -(A5(2:Nx+1,2:Ny+1,2:Nz+1) - A5(2:Nx+1,2:Ny+1,1:Nz)) / dz + sigma * ...
        ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A5(2:Nx+1,2:Ny+1,2:Nz+1));
    resA6 = (A6(2:Nx+1,2:Ny+1,3:Nz+2) - A6(2:Nx+1,2:Ny+1,2:Nz+1)) / dz + sigma * ...
        ((1/6)*(A1(2:Nx+1,2:Ny+1,2:Nz+1) + A2(2:Nx+1,2:Ny+1,2:Nz+1) + A3(2:Nx+1,2:Ny+1,2:Nz+1) + A4(2:Nx+1,2:Ny+1,2:Nz+1)+A5(2:Nx+1,2:Ny+1,2:Nz+1)+A6(2:Nx+1,2:Ny+1,2:Nz+1))-A6(2:Nx+1,2:Ny+1,2:Nz+1));
    
    res=norm(resA1, 'fro')^2 + norm(resA2, 'fro')^2 + norm(resA3, 'fro')^2 + norm(resA4, 'fro')^2 + norm(resA5, 'fro')^2 + norm(resA6, 'fro')^2;


    xslice=0;
    yslice=0;
    zslice=0;
    slice(X,Y,Z,A1(2:end-1,2:end-1,2:end-1)+A2(2:end-1,2:end-1,2:end-1)+A3(2:end-1,2:end-1,2:end-1)+A4(2:end-1,2:end-1,2:end-1)+A5(2:end-1,2:end-1,2:end-1)+A6(2:end-1,2:end-1,2:end-1),...
        xslice,yslice,zslice,'nearest');
    colorbar;

    pause(.1);

    n_iter = n_iter + 1;
    iter_arr=[iter_arr n_iter];
    res_arr=[res_arr,res];
    variableName=num2str(sigma);
    display=['sigma=',variableName];
    %plot(iter_arr,res_arr,'o-','LineWidth',2,'DisplayName',display);
    %xlabel('Number of iterations');
    %ylabel('Residual');
    %title('Convergence of the method');
    %grid on;
    %legend('show');

    fprintf("Iter %d: Residual: %f, Averge number of Newton iterations: %f\n", n_iter, res, n_Newton_iter/n_nonlin_eqs);
end
end