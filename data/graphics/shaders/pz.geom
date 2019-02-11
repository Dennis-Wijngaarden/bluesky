#version 330
#define VCOUNT_PZ 36
#define NM2M 1853.184
#define PI 3.1415926535897932384626433832795
layout (points) in;
layout (line_strip, max_vertices = 37) out;

// Uniform block of global data
layout (std140) uniform global_data {
int wrap_dir;           // Wrap-around direction
float wrap_lon;         // Wrap-around longitude
float panlat;           // Map panning coordinates [deg]
float panlon;           // Map panning coordinates [deg]
float zoom;             // Screen zoom factor [-]
int screen_width;       // Screen width in pixels
int screen_height;      // Screen height in pixels
int vertex_scale_type;  // Vertex scale type
};

in GSData {
    vec2 vAR;
    float rpz;
    vec4 color;
} gs_in[];

out vec4 color_fs;

void main()
{
    vec2 point_coord;
    for (int i = 0; i < VCOUNT_PZ + 1; i++)
    {
        color_fs = gs_in[0].color;
        point_coord = vec2(gs_in[0].rpz * cos(i * (2*PI) / VCOUNT_PZ), gs_in[0].rpz * sin(i * (2*PI) / VCOUNT_PZ));
        gl_Position = gl_in[0].gl_Position + vec4(gs_in[0].vAR * point_coord / NM2M, 0.0, 0.0) * zoom * 30.;
        EmitVertex();
    }
    EndPrimitive();
}