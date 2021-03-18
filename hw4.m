clear all; close all; clc

[images_train, labels_train] = mnist_parse('train-images-idx3-ubyte', 'train-labels-idx1-ubyte');
[images_test, labels_test] = mnist_parse('t10k-images-idx3-ubyte', 't10k-labels-idx1-ubyte');

images_train = im2double(images_train);
size_train = size(images_train);
images_train = reshape(images_train,size_train(1)*size_train(2),size_train(3));
images_train = images_train-repmat(mean(images_train,2),1,size_train(3));
[U,S,V] = svd(images_train,'econ');

sig = diag(S);
lambda = sig.^2;
threshold = 0.9;
energy = 0;
r = 0;
while energy <= threshold
    r = r+1;
    energy = energy+lambda(r)/sum(lambda);
end
figure(1)
subplot(2,1,1)
plot(lambda,'ko','Linewidth',1)
ylabel('\lambda');
xlabel('Measurements');
title("Singular Value Spectrum");
subplot(2,1,2)
semilogy(lambda,'ko','Linewidth',1)
ylabel('\lambda (log scale)');
xlabel('Measurements');
title("Singular Value Spectrum");

figure(2)
projection = U(:,[1,2,3])'*images_train;
projection0 = projection(:,find(labels_train == 0));
projection1 = projection(:,find(labels_train == 1));
projection2 = projection(:,find(labels_train == 2));
projection3 = projection(:,find(labels_train == 3));
projection4 = projection(:,find(labels_train == 4));
projection5 = projection(:,find(labels_train == 5));
projection6 = projection(:,find(labels_train == 6));
projection7 = projection(:,find(labels_train == 7));
projection8 = projection(:,find(labels_train == 8));
projection9 = projection(:,find(labels_train == 9));
plot3(projection0(1,:),projection0(2,:),projection0(3,:),'o');hold on
plot3(projection1(1,:),projection1(2,:),projection1(3,:),'o');hold on
plot3(projection2(1,:),projection2(2,:),projection2(3,:),'o');hold on
plot3(projection3(1,:),projection3(2,:),projection3(3,:),'o');hold on
plot3(projection4(1,:),projection4(2,:),projection4(3,:),'o');hold on
plot3(projection5(1,:),projection5(2,:),projection5(3,:),'o');hold on
plot3(projection6(1,:),projection6(2,:),projection6(3,:),'o');hold on
plot3(projection7(1,:),projection7(2,:),projection7(3,:),'o');hold on
plot3(projection8(1,:),projection8(2,:),projection8(3,:),'o');hold on
plot3(projection9(1,:),projection9(2,:),projection9(3,:),'o');
legend('0','1','2','3','4','5','6','7','8','9');
xlabel('Mode1');
ylabel('Mode2');
zlabel('Mode3');
title('Projection onto three selected V-modes: Column1,2,3');

% LDA-2 digits
feature = 12;
images = S*V';
i1 = find(labels_train == 1);
i2 = find(labels_train == 2);
dig1 = images(1:feature,i1);
dig2 = images(1:feature,i2);
num1 = length(i1);
num2 = length(i2);
% num1 = size(i1',2);
% num2 = size(i2',2);
m1 = mean(dig1,2);
m2 = mean(dig2,2);
Sw = 0; % within class variances
for k = 1:num1
    Sw = Sw + (dig1(:,k) - m1)*(dig1(:,k) - m1)';
end
for k = 1:num2
   Sw =  Sw + (dig2(:,k) - m2)*(dig2(:,k) - m2)';
end
Sb = (m1-m2)*(m1-m2)'; % between class

[V2, D] = eig(Sb,Sw); % linear disciminant analysis
[lambda2, ind] = max(abs(diag(D)));
w2 = V2(:,ind);
w2 = w2/norm(w2,2);

vdig1 = w2'*dig1;
vdig2 = w2'*dig2;

if mean(vdig1) > mean(vdig2)
    w2 = -w2;
    vdig1 = -vdig1;
    vdig2 = -vdig2;
end

figure(3)
plot(vdig1,0,'ob','Linewidth',2)
hold on
plot(vdig2,1,'dr','Linewidth',2)
ylabel('Specified value for separation');
xlabel('Sorted data');
title('Projection of the picked two digits');

sortdig1 = sort(vdig1);
sortdig2 = sort(vdig2);

t1 = length(sortdig1);
t2 = 1;
while sortdig1(t1) > sortdig2(t2)
    t1 = t1 - 1;
    t2 = t2 + 1;
end
threshold = (sortdig1(t1) + sortdig2(t2))/2;

figure(4)
subplot(1,2,1)
histogram(sortdig1,30); hold on, plot([threshold threshold],[0 1600],'r')
set(gca,'Ylim',[0 1600],'Fontsize',14);
title('Digit 1: ones');
subplot(1,2,2)
histogram(sortdig2,30); hold on, plot([threshold threshold],[0 700],'r')
set(gca,'Ylim',[0 700],'Fontsize',14)
title('Digit 2: twos')

% LDA-3 digits
i3 = find(labels_train == 1);
i4 = find(labels_train == 2);
i5 = find(labels_train == 3);
dig3 = images(1:feature,i3);
dig4 = images(1:feature,i4);
dig5 = images(1:feature,i5);
num3 = length(i3);
num4 = length(i4);
num5 = length(i5);

m3 = mean(dig3,2);
m4 = mean(dig4,2);
m5 = mean(dig5,2);
mj = [m3,m4,m5];
Sw = 0; % within class variances
for k = 1:num3
    Sw = Sw + (dig3(:,k) - m3)*(dig3(:,k) - m3)';
end
for k = 1:num4
   Sw =  Sw + (dig4(:,k) - m4)*(dig4(:,k) - m4)';
end
for k = 1:num5
   Sw =  Sw + (dig5(:,k) - m5)*(dig5(:,k) - m5)';
end

mn = (m3 + m4 + m5)/3;
Sb = (mj(:,1)-mn)*(mj(:,1)-mn)'+(mj(:,2)-mn)*(mj(:,2)-mn)'+(mj(:,3)-mn)*(mj(:,3)-mn)';

[V2, D] = eig(Sb,Sw); % linear disciminant analysis
[lambda3, ind] = max(abs(diag(D)));
w3 = V2(:,ind);
w3 = w3/norm(w3,2);

vdig3 = w3'*dig3;
vdig4 = w3'*dig4;
vdig5 = w3'*dig5;

if mean(vdig3) > mean(vdig4)
    w3 = -w3;
    vdig3 = -vdig3;
    vdig4 = -vdig4;
end

if mean(vdig4) > mean(vdig5) 
    w3 = -w3;
    vdig4 = -vdig4;
    vdig5 = -vdig5;
end

sortdig3 = sort(vdig3);
sortdig4 = sort(vdig4);
sortdig5 = sort(vdig5);

t3 = length(sortdig3);
t4 = 1;
while sortdig3(t3) > sortdig4(t4)
    t3 = t3 - 1;
    t4 = t4 + 1;
end
threshold1 = (sortdig3(t3) + sortdig4(t4))/2;

t3 = length(sortdig3);
t5=1;
while sortdig3(t3) > sortdig5(t5)
    t3 = t3 - 1;
    t5 = t5 + 1;
end
threshold2 = (sortdig3(t3) + sortdig5(t5))/2;
t4 = length(sortdig4);
t5=1;
while sortdig4(t4) > sortdig5(t5)
    t4 = t4 - 1;
    t5 = t5 + 1;
end
threshold3 = (sortdig4(t4) + sortdig5(t5))/2;
thresholdall = (threshold1+threshold2+threshold3)/3;

figure(5)
subplot(1,3,1)
histogram(sortdig3,30); hold on, plot([thresholdall thresholdall],[0,1100],'r')
set(gca,'Ylim',[0 1100],'Fontsize',14)
title('Digit 1: ones')
subplot(1,3,2)
histogram(sortdig4,30); hold on, plot([thresholdall thresholdall],[0,800],'r')
set(gca,'Ylim',[0 800],'Fontsize',14)
title('Digit 2: twos')
subplot(1,3,3)
histogram(sortdig5,30); hold on, plot([thresholdall thresholdall],[0,700],'r')
set(gca,'Ylim',[0 700],'Fontsize',14)
title('Digit 3: threes')

% find most difficult and easy
images_test = im2double(images_test);
size_test = size(images_test);
images_test = reshape(images_test,size_test(1)*size_test(2),size_test(3));
images_test = images_test-repmat(mean(images_test,2),1,size_test(3));
feature = 12;
nummax = 0;
nummin = 1;
max1 = 0;
max2 = 0;
min1 = 0;
min2 = 0;
for i = 1:9
    for j = (i + 1):10
        index1 = find(labels_train == i-1);
        index2 = find(labels_train == j-1);
        tdig1 = images_train(1:feature,index1);
        tdig2 = images_train(1:feature,index2);
        tnum1 = size(index1',2);
        tnum2 = size(index2',2);
        tm1 = mean(tdig1,2);
        tm2 = mean(tdig2,2);
        tSw = 0; % within class variances
        for k = 1:tnum1
            tSw = tSw + (tdig1(:,k) - tm1)*(tdig1(:,k) - tm1)';
        end
        for k = 1:tnum2
            tSw =  tSw + (tdig2(:,k) - tm2)*(tdig2(:,k) - tm2)';
        end
        tSb = (tm1-tm2)*(tm1-tm2)'; % between class

        [tV2, tD] = eig(tSb,tSw); % linear disciminant analysis
        [tlambda, tind] = max(abs(diag(tD)));
        tw = tV2(:, tind);
        tw = tw/norm(tw,2);

        tvdig1 = tw'*tdig1;
        tvdig2 = tw'*tdig2;
        
        if mean(tvdig1) > mean(tvdig2)
            tw = -tw;
            tvdig1 = -tvdig1;
            tvdig2 = -tvdig2;
        end

        tsortdig1 = sort(tvdig1);
        tsortdig2 = sort(tvdig2);

        tt1 = length(tsortdig1);
        tt2 = 1;
        while tsortdig1(tt1) > tsortdig2(tt2)
            tt1 = tt1 - 1;
            tt2 = tt2 + 1;
        end
        tthreshold = (tsortdig1(tt1) + tsortdig2(tt2))/2;
        
        TestNum = size(images_test,2);
        TestMat = U(:, 1:feature)'*images_test; % PCA projection
        pval = tw'*TestMat;
        i1test = find(labels_test == i-1);
        i2test = find(labels_test == j-1);
        ResVec1 = (pval(i1test) < threshold);
        err1 = abs(sum(ResVec1) - size(i1test',2));
        ResVec2 = (pval(i2test) > threshold);
        err2 = abs(sum(ResVec2) - size(i2test',2));        
        sucRate = 1 - (err1+err2)/TestNum;
        if sucRate > nummax
            nummax = sucRate;
            max1 = i-1;
            max2 = j-1;
        end
        if sucRate < nummin  
            nummin = sucRate;
            min1 = i-1;
            min2 = j-1;
        end     
    end
end

success_svm = 0;
Mdl = fitcecoc(images_train',labels_train);
test_labels = predict(Mdl,images_test');
for j = 1:length(test_labels)
    if (test_labels(j) - labels_test(j)) == 0
        success_svm = success_svm + 1;
    end
end
sucRate_svm = success_svm/length(test_labels);

tree = fitctree(images_train',labels_train, 'MaxNumSplits',10,'CrossVal','on');
view(tree.Trained{1},'Mode','graph');
classError = kfoldLoss(tree, 'mode', 'individual');
[~, k] = min(classError);
testpredict = predict(tree.Trained{k}, images_test');
errortree = immse(testpredict, labels_test);
diff = testpredict - labels_test;
success = find(diff == 0);
sucRate_tree = length(success) / length(labels_test);




