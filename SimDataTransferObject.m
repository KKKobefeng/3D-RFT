classdef SimDataTransferObject
    %SimDataTransferObject

    properties
        points;
        c_inc;
        v_norm_vec; 
        F; 
        f;
        forces_x; 
        forces_y; 
        forces_z; 
        forces; 
        T;
        torque_x; 
        torque_y; 
        torque_z;
        alpha_gen;
        alpha_gen_n; 
        alpha_gen_t; 
        alpha;
        TRG; 
        TRG_visual;
    end

    methods
        function obj = SimDataTransferObject(points, c_inc, v_norm_vec, F, f, forces_x, ...
                forces_y, forces_z, forces, T, torque_x, torque_y, torque_z, ...
                alpha_gen, alpha_gen_n, alpha_gen_t, alpha, TRG, TRG_visual)
            obj.points = points;
            obj.c_inc = c_inc;
            obj.v_norm_vec = v_norm_vec;
            obj.F = F;
            obj.f = f;
            obj.forces_x = forces_x;
            obj.forces_y = forces_y;
            obj.forces_z = forces_z;
            obj.forces = forces;
            obj.T = T;
            obj.torque_x = torque_x;
            obj.torque_y = torque_y;
            obj.torque_z = torque_z;
            obj.alpha_gen = alpha_gen;
            obj.alpha_gen_n = alpha_gen_n;
            obj.alpha_gen_t = alpha_gen_t;
            obj.alpha = alpha;
            obj.TRG = TRG;
            obj.TRG_visual = TRG_visual;
        end
    end
end