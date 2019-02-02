#version 330
in vec2 ssd_coord;
in vec4 color_fs;
in vec4 asasown_frag;
out vec4 color_out;

// Vlimits is [Vmin^2, Vmax^2, Vmax]
uniform vec3 Vlimits;

void main()
{
    float Vsquared    = dot(ssd_coord, ssd_coord);
    if (Vsquared < asasown_frag[1] || Vsquared > asasown_frag[2]) {
        discard;
    } else {
        //color_out = vec4(color_fs, smoothstep(asasown_frag[1], asasown_frag[2], Vsquared));
        color_out = color_fs;
    }
}
