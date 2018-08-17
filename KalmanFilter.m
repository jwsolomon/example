%Kalman Filtering

measurements = load('measurements.txt');
trajectory = load('trajectory.txt');

%Particle filter

Ns = 1000;
Nt = 100;

P0 = 2;
x0 = 4;
Q = 1;
R = 0.25;
N = 100;

%Script

P = P0;
x_start = x0;
[ws, logweights] = resetWeights();
particles_k = resetParticles(x_start, P);
noise = resetNoise();
Neff = Ns;

P_array = zeros(N, 1);
x_array = zeros(N, 1);
P_array(1) = P;
x_array(1,:) = x_start;
for i = 1:N-1
    new_particles = propagateParticles(particles_k, noise, i);
    [logweights, ws] = updateWeights(i, new_particles, logweights);
    [x_new, P_new] = computeEstimate(ws, new_particles);
    Neff = effectiveParticles(ws);
    
    P_array(i) = P_new;
    x_array(i) = x_new;
    
    if Neff < Nt
        particles_k = resetParticles(x_new, P_new);
        noise = resetNoise();
        [ws, logweights] = resetWeights();
    else
        particles_k = new_particles;
        noise = resetNoise();
    end
end
        
        


%function definitions


function exup = expectedUpdate(x,k)
    exup = 2*atan(x) + 0.5*cos((pi/3)*k);
end

function exmes = expectedMeasurement(x,k)
    exmes = x + x^2 + x^3;
end

function [ws,logweights] = resetWeights()
    Ns = 1000;
    ws = ones(Ns,1)*(1./Ns);
    logweights = log(ws);
end

function particles_k = resetParticles(x,P)
    Ns = 1000;
    particles_k = mvnrnd(x,P,Ns);
end

function noise_k = resetNoise()
    Q = 1;
    Ns = 1000;
    noise_k = mvnrnd(zeros(1),Q,Ns);
end

function new_particles = propagateParticles(particles,noise,k)
    Ns = 1000;
    new_particles = zeros(Ns,1);
    for i = 1:Ns
        new_particles(i,:) = expectedUpdate(particles(i,1),k) + noise(i,1);
    end
end

function [logweights,new_weights] = updateWeights(index, new_particles, logweights)
    measurements = load('measurements.txt');
    R = 0.25;
    measurement = measurements(index);
    Ns = 1000;
    for i = 1:length(Ns)
        logweights(i) = (-.5*(measurement - expectedMeasurement(new_particles(i)))*(1./R)*(measurement - expectedMeasurement(new_particles(i)))+ logweights(i));
        maxlog = max(logweights);
        intermediate_weights = exp(logweights - maxlog);
        new_weights = intermediate_weights/(sum(intermediate_weights));
    end
end


function [x_new,P_new] = computeEstimate(weights, particles)
    Ns = 1000;
    x_new = zeros(1);
    P_new = zeros(1);
    for i = 1:length(Ns)
        x_new = weights(i)*particles(i) + x_new;
    end
    for i = 1:length(Ns)
        diff = particles(i) - x_new;
        P_new = weights(i)*dot(diff', diff);
    end
end


function efpart = effectiveParticles(weights)
    efpart = 1/(sum(weights.*weights));
end