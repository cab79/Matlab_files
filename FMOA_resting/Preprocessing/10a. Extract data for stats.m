clear all

subjects = {'H1', 'H2', 'H3', 'H4', 'H5', 'H6', 'H7', 'H8', 'H9', 'H10', 'H11', 'H12', 'H13' 'H14', 'H15', 'H16', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7', 'F8', 'F9', 'F10', 'F11', 'F12', 'F13', 'F14', 'F15', 'F16', 'OA1', 'OA2', 'OA4', 'OA5', 'OA6', 'OA7', 'OA8', 'OA9', 'OA10', 'OA11', 'OA12', 'OA13', 'OA14', 'OA15', 'OA16', 'OA17'};

load gamma1;
load gamma2;
load gamma3;
load beta;
load alpha1;
load alpha2;
load alpha3;
load theta;
load delta;

gamma1_open = gamma1(:,:,1);
save gamma1_open.mat gamma1_open
gamma1_closed = gamma1(:,:,2);
save gamma1_closed.mat gamma1_closed

gamma2_open = gamma2(:,:,1);
save gamma2_open.mat gamma2_open
gamma2_closed = gamma2(:,:,2);
save gamma2_closed.mat gamma2_closed

gamma3_open = gamma3(:,:,1);
save gamma3_open.mat gamma3_open
gamma3_closed = gamma3(:,:,2);
save gamma3_closed.mat gamma3_closed

beta_open = beta(:,:,1);
save beta_open.mat beta_open
beta_closed = beta(:,:,2);
save beta_closed.mat beta_closed

alpha1_open = alpha1(:,:,1);
save alpha1_open.mat alpha1_open
alpha1_closed = alpha1(:,:,2);
save alpha1_closed.mat alpha1_closed

alpha2_open = alpha2(:,:,1);
save alpha2_open.mat alpha2_open
alpha2_closed = alpha2(:,:,2);
save alpha2_closed.mat alpha2_closed

alpha3_open = alpha3(:,:,1);
save alpha3_open.mat alpha3_open
alpha3_closed = alpha3(:,:,2);
save alpha3_closed.mat alpha3_closed

theta_open = theta(:,:,1);
save theta_open.mat theta_open
theta_closed = theta(:,:,2);
save theta_closed.mat theta_closed

delta_open = delta(:,:,1);
save delta_open.mat delta_open
delta_closed = delta(:,:,2);
save delta_closed.mat delta_closed

