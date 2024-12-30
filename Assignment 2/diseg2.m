function diseg2(modes,scale_factor,incidenze,l,gamma,posiz,idb,xy,margin_rel)
if ~exist('margin_rel', 'var')
    margin_rel = 0.1;
end

% Checking consistency input data
[n_el,~]=size(incidenze);

if length(posiz) ~= n_el
    error('Error: number of nodes in posiz matrix differs from number of elements')
end

n_gdl=size(modes,1);
n_modes=size(modes,2);


hold on
colors = lines(n_modes);
LH = gobjects(1, n_modes + 1); % Preallocate line handles array
L = cell(1, n_modes + 1); % Preallocate legend labels cell array

% Plot the reference structure (undeformed)
grey = [0.3, 0.3, 0.3];
for k = 1:n_el
    % Computing the axial (u) and transversal (w) displacements by means of shape functions
    csi = l(k) * (0:0.05:1);
    lambda = [cos(gamma(k)) sin(gamma(k)) 0.
             -sin(gamma(k)) cos(gamma(k)) 0.
              0.           0.           1.];
    undef = lambda(1:2, 1:2)' * [csi; zeros(1, length(csi))];
    plot(undef(1, :) + posiz(k, 1), undef(2, :) + posiz(k, 2), '--', 'LineWidth', 1.5, 'HandleVisibility', 'off', 'Color', grey);
end
plot(xy(:, 1), xy(:, 2), '.', 'MarkerSize', 15, 'HandleVisibility', 'off', 'Color', grey);
LH(1) = plot(nan, nan, '.--', 'LineWidth', 1.5, 'MarkerSize', 15, 'HandleVisibility', 'off', 'Color', grey);
L{1} = 'Reference structure';
    

% Looping on the modes
for mode_idx = 1:n_modes
    mode = modes(:,mode_idx);
    color = colors(mode_idx, :);

    % looping on the finite elements
    for k=1:n_el
        % building the nodal displacements vector of each element in the global reference frame
        xkG=zeros(6,1);
        for iri=1:6
            if incidenze(k,iri) <= n_gdl
                xkG(iri,1)=mode(incidenze(k,iri));
            else
                xkG(iri,1)=0.;
            end
        end
        % Applying scale factor
        xkG=scale_factor*xkG;
        % Global to Local reference frame rotation
        lambda = [cos(gamma(k)) sin(gamma(k)) 0. 
                 -sin(gamma(k)) cos(gamma(k)) 0.
                  0.            0.            1.];
        Lambda = [lambda     zeros(3)
                  zeros(3)   lambda];
        xkL = Lambda * xkG;

        % Computing the axial (u) and transversal (w) displacements by means of shape functions
        csi=l(k)*(0:0.05:1);
        fu=zeros(6,length(csi));
        fu(1,:)=1-csi/l(k);
        fu(4,:)=csi/l(k);
        u=(fu'*xkL)';
        fw=zeros(6,length(csi));
        fw(2,:)=2*(csi/l(k)).^3-3*(csi/l(k)).^2+1;
        fw(3,:)=l(k)*((csi/l(k)).^3-2*(csi/l(k)).^2+csi/l(k));
        fw(5,:)=-2*(csi/l(k)).^3+3*(csi/l(k)).^2;
        fw(6,:)=l(k)*((csi/l(k)).^3-(csi/l(k)).^2);
        w=(fw'*xkL)';
        % Local to global transformation of the element's deformation
        xyG=lambda(1:2,1:2)'*[u+csi;w];
        % Plotting deformed element's shape
        plot(xyG(1,:)+posiz(k,1),xyG(2,:)+posiz(k,2),'Color', color, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    end

    % Looping the nodes
    n_nodi=size(idb,1);
    xkG=zeros(n_nodi,2);
    for k=1:n_nodi
        for ixy=1:2
            if idb(k,ixy) <= n_gdl
                xkG(k,ixy)=mode(idb(k,ixy));
            else
                xkG(k,ixy)=0.;
            end
        end
    end
    xkG=scale_factor*xkG;
    xyG=xkG+xy;
    plot(xyG(:,1),xyG(:,2),'o', 'Color', color,'LineWidth', 1.5, 'HandleVisibility', 'off');
    LH(mode_idx + 1) = plot(nan, nan, 'o-', 'Color', color, 'LineWidth', 1.5, 'HandleVisibility', 'off');
    L{mode_idx + 1} = ['Mode ', num2str(mode_idx)];
end

grid on
box on
axis equal

% Calculate the limits for the plot with a slight margin
x_min = min(xy(:,1));
x_max = max(xy(:,1));
y_min = min(xy(:,2));
y_max = max(xy(:,2));
margin = margin_rel * max(x_max - x_min, y_max - y_min);
xlim([x_min - margin, x_max + margin]);
ylim([y_min - margin, y_max + margin]);

xlabel('x [m]')
ylabel('y [m]')

if n_modes <= 1
    L{2} = 'Mode shape';
end
legend(LH, L, 'Location', 'Best')