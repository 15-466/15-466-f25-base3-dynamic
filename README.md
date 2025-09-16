# Dynamic Meshes with Scene

This repository demonstrates how to use dynamic geometry with `Scene`.

It introduces a new struct `DynamicMeshBuffer` which manages a vertex buffer.

The recipe here is:

1. Add a `DynamicMeshBuffer` as a member variable of your mode.
2. In the mode's construtor, build a vertex array object that maps the attributes held in the buffer to programs you are going to use to draw them.
3. At any time, create `Scene::Drawable`s that reference the `DynamicMeshBuffer` via the vao created in step 2.
4. At any time, update the contents of the `DynamicMeshBuffer` with its `set()` member function. (Remember to update the `count` of any drawables that reference it as well!)
5. In the mode's destructor, remember to delete the vertex array object.

Note: there isn't much code here, and it maps pretty directly to an OpenGL vertex buffer -- `set` is basically `glBufferData`, and the constructor and destructor are `glGenBuffers` / `glDeleteBuffers`. The only slightly complicated function is the `make_vao_for_program` member, because that looks up attribute locations.

Note2: you could potentially save (several!) copies by using `glMapBuffer` and generating geometry directly into mapped buffer member instead of `glBufferData` which must copy out of the passed pointer.
