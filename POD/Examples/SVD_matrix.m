%% Compute SVD
close all
clear all
clc
matrix
[U,S,V] = svd(Z);
figure
subplot(2,2,1)
surf(X,T,Z, 'Linestyle','none')
axis([0,1,0,2,0.4,2.1])
xlabel('x'), ylabel('y'), zlabel('z'), title('Actual surface')


for k = 1:3
    ZZ = U(:,1:k) * S(1:k,1:k) * V(:,1:k)';
    subplot(2,2,k+1)
    surf(X,T,ZZ, 'Linestyle','none')
    axis([0,1,0,2,0.4,2.1])
    xlabel('x'), ylabel('y'), zlabel('z'), title(['Rank ' num2str(k) ' approximation'])
end
figure
subplot(1,2,1)
S = diag(S);
semilogy(S, 'o')
axis tight
xlabel('number'), ylabel('singular value'), title(['Singular Values of Z'])


subplot(1,2,2)
V1 = V(:,1:3);
C1 = V1'* Z';

plot(t, C1(1,:), '--r', t, C1(2,:), ':b', t, C1(3,:), '-.g')
xlabel('t'), ylabel('Modal coordinates'), title(['Modal contributions'])
legend('1','2','3')


C = V'* Z';
figure
mesh(C)
figure
for i = 1 : length(S)
    plot(t,C(i,:))
    hold on
end

