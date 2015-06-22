#include "Scene_polylines_item.h"

#include <CGAL/bounding_box.h>
#include <CGAL/gl.h>
#include <CGAL/glu.h>
#include <QMenu>
#include <QAction>

#include <QInputDialog>
namespace {
void CGALglcolor(QColor c, int dv = 0)
{
    if ( 0 != dv )
    {
        // workaround for Qt-4.2.
#if QT_VERSION < 0x040300
#  define darker dark
#endif
        c = c.darker(dv);
#undef darker
    }
    //::glColor4d(c.red()/255.0, c.green()/255.0, c.blue()/255.0, c.alpha()/255.0);
}
}
struct light_info
{
    //position
    GLfloat position[4];

    //ambient
    GLfloat ambient[4];

    //diffuse
    GLfloat diffuse[4];

    //specular
    GLfloat specular[4];
};
typedef Scene_polylines_item::K K;
typedef K::Point_3 Point_3;
//Fill the VBO with coordinates of the vertices composing a sphere
void Scene_polylines_item::create_Sphere(double R)
{
    float T, P;
    float x[4],y[4],z[4];


    //Top of the sphere
    for(int t=0; t<360; t+=sectors)
    {

        positions_spheres.push_back(0);
        positions_spheres.push_back(0);
        positions_spheres.push_back(R);
        positions_spheres.push_back(1.0);


        normals_spheres.push_back(0);
        normals_spheres.push_back(0);
        normals_spheres.push_back(1);



        P = rings*M_PI/180.0;
        T = t*M_PI/180.0;
        x[1] = sin(P) * cos(T) ;
        y[1] = sin(P) * sin(T) ;
        z[1] = cos(P);
        positions_spheres.push_back(R * x[1]);
        positions_spheres.push_back(R * y[1]);
        positions_spheres.push_back(R * z[1]);
        positions_spheres.push_back(1.0);

        normals_spheres.push_back(x[1]);
        normals_spheres.push_back(y[1]);
        normals_spheres.push_back(z[1]);

        //
        P = rings*M_PI/180.0;
        T = (t+sectors)*M_PI/180.0;
        x[2] = sin(P) * cos(T) ;
        y[2] = sin(P) * sin(T) ;
        z[2] = cos(P);
        positions_spheres.push_back(R * x[2]);
        positions_spheres.push_back(R * y[2]);
        positions_spheres.push_back(R * z[2]);
        positions_spheres.push_back(1.0);

        normals_spheres.push_back(x[2]);
        normals_spheres.push_back(y[2]);
        normals_spheres.push_back(z[2]);

    }

    //Body of the sphere
    for (int p=rings; p<180-rings; p+=rings)
        for(int t=0; t<360; t+=sectors)
        {
            //A
            P = p*M_PI/180.0;
            T = t*M_PI/180.0;
            x[0] = sin(P) * cos(T) ;
            y[0] = sin(P) * sin(T) ;
            z[0] = cos(P);

            positions_spheres.push_back(R * x[0]);
            positions_spheres.push_back(R * y[0]);
            positions_spheres.push_back(R * z[0]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[0]);
            normals_spheres.push_back(y[0]);
            normals_spheres.push_back(z[0]);

            //B
            P = (p+rings)*M_PI/180.0;
            T = t*M_PI/180.0;
            x[1] = sin(P) * cos(T) ;
            y[1] = sin(P) * sin(T) ;
            z[1] = cos(P);
            positions_spheres.push_back(R * x[1]);
            positions_spheres.push_back(R * y[1]);
            positions_spheres.push_back(R * z[1]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[1]);
            normals_spheres.push_back(y[1]);
            normals_spheres.push_back(z[1]);

            //C
            P = p*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[2] = sin(P) * cos(T) ;
            y[2] = sin(P) * sin(T) ;
            z[2] = cos(P);
            positions_spheres.push_back(R * x[2]);
            positions_spheres.push_back(R * y[2]);
            positions_spheres.push_back(R * z[2]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[2]);
            normals_spheres.push_back(y[2]);
            normals_spheres.push_back(z[2]);
            //D
            P = (p+rings)*M_PI/180.0;
            T = (t+sectors)*M_PI/180.0;
            x[3] = sin(P) * cos(T) ;
            y[3] = sin(P) * sin(T) ;
            z[3] = cos(P);
            positions_spheres.push_back(R * x[3]);
            positions_spheres.push_back(R * y[3]);
            positions_spheres.push_back(R * z[3]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[3]);
            normals_spheres.push_back(y[3]);
            normals_spheres.push_back(z[3]);



            positions_spheres.push_back(R * x[1]);
            positions_spheres.push_back(R * y[1]);
            positions_spheres.push_back(R * z[1]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[1]);
            normals_spheres.push_back(y[1]);
            normals_spheres.push_back(z[1]);

            positions_spheres.push_back(R * x[2]);
            positions_spheres.push_back(R * y[2]);
            positions_spheres.push_back(R * z[2]);
            positions_spheres.push_back(1.0);

            normals_spheres.push_back(x[2]);
            normals_spheres.push_back(y[2]);
            normals_spheres.push_back(z[2]);

        }
    //Bottom of the sphere
    for(int t=0; t<360; t+=sectors)
    {


        positions_spheres.push_back(0);
        positions_spheres.push_back(0);
        positions_spheres.push_back(-R);
        positions_spheres.push_back(1.0);

        normals_spheres.push_back(0);
        normals_spheres.push_back(0);
        normals_spheres.push_back(-1);


        P = (180-rings)*M_PI/180.0;
        T = t*M_PI/180.0;
        x[1] = sin(P) * cos(T) ;
        y[1] = sin(P) * sin(T) ;
        z[1] = cos(P);
        positions_spheres.push_back(R * x[1]);
        positions_spheres.push_back(R * y[1]);
        positions_spheres.push_back(R * z[1]);
        positions_spheres.push_back(1.0);

        normals_spheres.push_back(x[1]);
        normals_spheres.push_back(y[1]);
        normals_spheres.push_back(z[1]);


        P = (180-rings)*M_PI/180.0;
        T = (t+sectors)*M_PI/180.0;
        x[2] = sin(P) * cos(T) ;
        y[2] = sin(P) * sin(T) ;
        z[2] = cos(P);
        positions_spheres.push_back(R * x[2]);
        positions_spheres.push_back(R * y[2]);
        positions_spheres.push_back(R * z[2]);
        positions_spheres.push_back(1.0);

        normals_spheres.push_back(x[2]);
        normals_spheres.push_back(y[2]);
        normals_spheres.push_back(z[2]);

    }
}

class Scene_polylines_item_private {
public:
    typedef Scene_polylines_item::K K;
    typedef K::Point_3 Point_3;

    Scene_polylines_item_private() :
        draw_extremities(false),
        spheres_drawn_radius(0),
        sphere_display_list(0)/*,
        quadric(0)*/
    {}

    ~Scene_polylines_item_private()
    {
      //  if(quadric != 0)
      //      gluDeleteQuadric(quadric);
      //  if(sphere_display_list  != 0)
      //      glDeleteLists(sphere_display_list, 1);
    }

    void draw_sphere(const K::Point_3&, double) const;
    void draw_spheres(const Scene_polylines_item*) const;

    bool draw_extremities;
    double spheres_drawn_radius;
private:
    mutable GLuint sphere_display_list;
   // mutable GLUquadric* quadric;
};

void
Scene_polylines_item::initialize_buffers(Viewer_interface *viewer = 0) const
{
//vao for the lines
    {
        program = getShaderProgram(PROGRAM_WITHOUT_LIGHT, viewer);
        program->bind();

        vaos[0]->bind();
        buffers[0].bind();
        buffers[0].allocate(positions_lines.data(), positions_lines.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[0].release();
        vaos[0]->release();
        program->release();
    }
   //vao for the spheres
    {
        program = getShaderProgram(PROGRAM_INSTANCED, viewer);
        program->bind();

        vaos[1]->bind();
        buffers[1].bind();
        buffers[1].allocate(positions_spheres.data(), positions_spheres.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[1].release();

        buffers[2].bind();
        buffers[2].allocate(normals_spheres.data(), normals_spheres.size()*sizeof(float));
        program->enableAttributeArray("normals");
        program->setAttributeBuffer("normals",GL_FLOAT,0,3);
        buffers[2].release();

        buffers[3].bind();
        buffers[3].allocate(color_spheres.data(), color_spheres.size()*sizeof(float));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[3].release();

        buffers[4].bind();
        buffers[4].allocate(positions_center.data(), positions_center.size()*sizeof(float));
        program->enableAttributeArray("center");
        program->setAttributeBuffer("center",GL_FLOAT,0,3);
        buffers[4].release();

        qFunc.glVertexAttribDivisor(program->attributeLocation("center"), 1);
        qFunc.glVertexAttribDivisor(program->attributeLocation("colors"), 1);
        vaos[1]->release();

        program->release();
    }

//vao for the wired spheres
    {
        program = getShaderProgram(PROGRAM_INSTANCED_WIRE, viewer);
        program->bind();

        vaos[2]->bind();
        buffers[5].bind();
        buffers[5].allocate(positions_wire_spheres.data(), positions_wire_spheres.size()*sizeof(float));
        program->enableAttributeArray("vertex");
        program->setAttributeBuffer("vertex",GL_FLOAT,0,4);
        buffers[5].release();
        QColor temp = this->color();
        program->setAttributeValue("colors", temp);

        program->setAttributeValue("normals",QVector3D(0.0,0.0,0.0));

        buffers[6].bind();
        buffers[6].allocate(color_spheres.data(), color_spheres.size()*sizeof(float));
        program->enableAttributeArray("colors");
        program->setAttributeBuffer("colors",GL_FLOAT,0,3);
        buffers[6].release();

        buffers[7].bind();
        buffers[7].allocate(positions_center.data(), positions_center.size()*sizeof(float));
        program->enableAttributeArray("center");
        program->setAttributeBuffer("center",GL_FLOAT,0,3);
        buffers[7].release();

        qFunc.glVertexAttribDivisor(program->attributeLocation("center"), 1);
        vaos[2]->release();
        program->release();
    }

   are_buffers_filled = true;

}
void
Scene_polylines_item::compute_elements()
{
    positions_spheres.clear();
    positions_wire_spheres.clear();
    positions_lines.clear();
    color_spheres.clear();
    normals_spheres.clear();
    positions_center.clear();
    nbSpheres = 0;

    //Fills the VBO with the lines
    for(std::list<std::vector<Point_3> >::const_iterator it = polylines.begin();
        it != polylines.end();
        ++it){
        if(it->empty()) continue;
        for(size_t i = 0, end = it->size()-1;
            i < end; ++i)
        {
            const Point_3& a = (*it)[i];
            const Point_3& b = (*it)[i+1];
            positions_lines.push_back(a.x());
            positions_lines.push_back(a.y());
            positions_lines.push_back(a.z());
            positions_lines.push_back(1.0);

            positions_lines.push_back(b.x());
            positions_lines.push_back(b.y());
            positions_lines.push_back(b.z());
            positions_lines.push_back(1.0);

        }

    }
    //Fills the VBO with the spheres
    if(d->draw_extremities)
    {

        // FIRST, count the number of incident cycles and polylines
        // for all extremities.
        typedef std::map<Point_3, int> Point_to_int_map;
        typedef Point_to_int_map::iterator iterator;
        Point_to_int_map corner_polyline_nb;

        { // scope to fill corner_polyline_nb'
            Point_to_int_map corner_cycles_nb;

            for(std::list<std::vector<Point_3> >::const_iterator
                it = this->polylines.begin(),
                end = this->polylines.end();
                it != end; ++it)
            {
                const K::Point_3& a = *it->begin();
                const K::Point_3& b = *it->rbegin();
                if(a == b) {
                    if ( it->size()>1 )
                        ++corner_cycles_nb[a];
                    else
                        ++corner_polyline_nb[a];
                }
                else {
                    ++corner_polyline_nb[a];
                    ++corner_polyline_nb[b];
                }
            }
            // THEN, ignore points that are incident to one cycle only.
            for(iterator
                c_it = corner_cycles_nb.begin(),
                end = corner_cycles_nb.end();
                c_it != end; ++c_it)
            {
                const Point_3& a = c_it->first;

                iterator p_it = corner_polyline_nb.find(a);

                // If the point 'a'=c_it->first has only incident cycles...
                if(p_it == corner_polyline_nb.end()) {
                    // ...then count it as a corner only if it has two incident cycles
                    // or more.
                    if(c_it->second > 1) {
                        corner_polyline_nb[a] = c_it->second;
                    }
                } else {
                    // else add the number of cycles.
                    p_it->second += c_it->second;
                }
            }
        }
        // At this point, 'corner_polyline_nb' gives the multiplicity of all
        // corners.
        //Finds the centers of the spheres and their color
        for(iterator
            p_it = corner_polyline_nb.begin(),
            end = corner_polyline_nb.end();
            p_it != end; ++p_it)
        {
            nbSpheres++;
            const K::Point_3& centre = p_it->first;
            positions_center.push_back(centre.x());
            positions_center.push_back(centre.y());
            positions_center.push_back(centre.z());

            float colors[3];
            switch(p_it->second) {
            case 1:
                colors[0] = 0.0; // black
                colors[1] = 0.0;
                colors[2] = 0.0;
                break;
            case 2:
                colors[0] = 0.0; // green
                colors[1] = 0.8;
                colors[2] = 0.0;
                break;
            case 3:
                colors[0] = 0.0; // blue
                colors[1] = 0.0;
                colors[2] = 0.8;
                break;
            case 4:
                colors[0] = 0.8; //red
                colors[1] = 0.0;
                colors[2] = 0.0;
                break;
            default:
                colors[0] = 0.8; //fuschia
                colors[1] = 0.0;
                colors[2] = 0.8;
            }

            color_spheres.push_back(colors[0]);
            color_spheres.push_back(colors[1]);
            color_spheres.push_back(colors[2]);
        }
        create_Sphere(d->spheres_drawn_radius);

        //Convert the triangle coordinates to lines coordinates for the
        //Wiremode in the spheres
        for(int i=0; i< positions_spheres.size(); i=i)
        {
            //draw triangles
            if(i< (360/sectors)*12)
            {
                //AB
                for(int j=i; j<i+8; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //BC
                for(int j=i+4; j<i+12; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //CA
                for(int j=i+8; j<i+16; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j%12]);
                }
                i+=12;
            }
            //draw quads
            else if((360/sectors) * 3 * 4 < i < positions_spheres.size() - (360/sectors) * 3 * 4)
            {
                //AB
                for(int j=i; j<i+8; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //BD
                for(int j=i+4; j<i+8; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                for(int j=i+12; j<i+16; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //DC
                for(int j=i+12; j<i+16; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                for(int j=i+8; j<i+12; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //CA
                for(int j=i+8; j<i+16; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j%12]);
                }
                i+=24;
            }
            //draw triangles
            else
            {
                //AB
                for(int j=i; j<i+8; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //BC
                for(int j=i+4; j<i+12; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j]);
                }
                //CA
                for(int j=i+8; j<i+16; j++)
                {
                    positions_wire_spheres.push_back(positions_spheres[j%12]);
                }
                i+=12;
            }

        }
    }


}


Scene_polylines_item::Scene_polylines_item() 
    : d(new Scene_polylines_item_private()),positions_lines(0), positions_spheres(0),
      normals_spheres(0), positions_center(0),color_spheres(0), positions_wire_spheres(0),nbSpheres(0),
      rings(18), sectors(36), Scene_item(8,3)
{
    qFunc.initializeOpenGLFunctions();
    changed();

}

Scene_polylines_item::~Scene_polylines_item()
{
    delete d;

}

bool
Scene_polylines_item::isEmpty() const {
    return polylines.empty();
}

Scene_interface::Bbox 
Scene_polylines_item::bbox() const {
    if(isEmpty())
        return Bbox();
    std::list<Point_3> boxes;
    for(std::list<std::vector<Point_3> >::const_iterator it = polylines.begin();
        it != polylines.end();
        ++it){
        if(it->begin() != it->end()) {
            Iso_cuboid_3 cub = CGAL::bounding_box(it->begin(), it->end());
            boxes.push_back((cub.min)());
            boxes.push_back((cub.max)());
        }
    }
    Iso_cuboid_3 bbox =
            boxes.begin() != boxes.end() ?
                CGAL::bounding_box(boxes.begin(), boxes.end()) :
                Iso_cuboid_3();

    return Bbox(bbox.xmin(),
                bbox.ymin(),
                bbox.zmin(),
                bbox.xmax(),
                bbox.ymax(),
                bbox.zmax());
}

Scene_polylines_item* 
Scene_polylines_item::clone() const {
    Scene_polylines_item* item = new Scene_polylines_item;
    item->polylines = polylines;
    QVariant metadata_variant = property("polylines metadata");
    if(metadata_variant.type() == QVariant::StringList)
    {
        item->setProperty("polylines metadata", metadata_variant);
    }
    return item;
}

QString
Scene_polylines_item::toolTip() const {
    QString s =
            tr("<p><b>%1</b> (mode: %2, color: %3)<br />"
               "<i>Polylines</i></p>"
               "<p>Number of polylines: %4</p>")
            .arg(this->name())
            .arg(this->renderingModeName())
            .arg(this->color().name())
            .arg(polylines.size());
    if(d->draw_extremities) {
        s += tr("<p>Legende of endpoints colors: <ul>"
                "<li>black: one incident polyline</li>"
                "<li>green: two incident polylines</li>"
                "<li>blue: three incident polylines</li>"
                "<li>red: four incident polylines</li>"
                "<li>fuchsia: five or more incident polylines</li>"
                "</ul></p>");
    }
    return s;
}

bool
Scene_polylines_item::supportsRenderingMode(RenderingMode m) const {
    return (m == Wireframe ||
            m == FlatPlusEdges ||
            m == Points);
}

// Shaded OpenGL drawing: only draw spheres
void
Scene_polylines_item::draw(Viewer_interface* viewer) const {

    if(!are_buffers_filled)
        initialize_buffers(viewer);
    if(d->draw_extremities)
    {
        vaos[1]->bind();
        program = getShaderProgram(PROGRAM_INSTANCED);
        attrib_buffers(viewer, PROGRAM_INSTANCED);
        program->bind();
        qFunc.glDrawArraysInstanced(GL_TRIANGLES, 0, positions_spheres.size()/4, nbSpheres);
        program->release();
        vaos[1]->release();
    }
}

// Wireframe OpenGL drawing
void 
Scene_polylines_item::draw_edges(Viewer_interface* viewer) const {

    if(!are_buffers_filled)
        initialize_buffers(viewer);

    vaos[0]->bind();
    attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    program->bind();
    QColor temp = this->color();
    program->setAttributeValue("colors", temp);
    qFunc.glDrawArrays(GL_LINES, 0, positions_lines.size()/4);
    program->release();
    vaos[0]->release();
    if(d->draw_extremities)
    {
        vaos[2]->bind();
        attrib_buffers(viewer, PROGRAM_INSTANCED_WIRE);
        program = getShaderProgram(PROGRAM_INSTANCED_WIRE);
        program->bind();
        qFunc.glDrawArraysInstanced(GL_LINES, 0, positions_wire_spheres.size()/4, nbSpheres);
        program->release();
        vaos[2]->release();
    }

}

void 
Scene_polylines_item::draw_points(Viewer_interface* viewer) const {
    if(!are_buffers_filled)
        initialize_buffers(viewer);

    vaos[0]->bind();
    attrib_buffers(viewer, PROGRAM_WITHOUT_LIGHT);
    program = getShaderProgram(PROGRAM_WITHOUT_LIGHT);
    program->bind();
    qFunc.glDrawArrays(GL_POINTS, 0, positions_lines.size()/4);
    // Clean-up
   vaos[0]->release();
   program->release();
}

QMenu* Scene_polylines_item::contextMenu() 
{
    const char* prop_name = "Menu modified by Scene_polylines_item.";

    QMenu* menu = Scene_item::contextMenu();

    // Use dynamic properties:
    // http://doc.trolltech.com/lastest/qobject.html#property
    bool menuChanged = menu->property(prop_name).toBool();

    if(!menuChanged) {
        menu->addSeparator();
        // TODO: add actions to display corners
        QAction* action = menu->addAction(tr("Display corners with radius..."));
        connect(action, SIGNAL(triggered()),
                this, SLOT(change_corner_radii()));

        QAction* actionSmoothPolylines =
                menu->addAction(tr("Smooth polylines"));
        actionSmoothPolylines->setObjectName("actionSmoothPolylines");
        connect(actionSmoothPolylines, SIGNAL(triggered()),this, SLOT(smooth()));
        menu->setProperty(prop_name, true);
    }
    return menu;
}

void Scene_polylines_item::changed()
{
    compute_elements();
    are_buffers_filled = false;


}

void Scene_polylines_item::change_corner_radii() {
    bool ok = true;
    double proposed_radius = d->spheres_drawn_radius;
    if(proposed_radius == 0) {
        Scene_interface::Bbox b = bbox();
        proposed_radius = (std::max)(b.xmax - b.xmin,
                                     proposed_radius);
        proposed_radius = (std::max)(b.ymax - b.ymin,
                                     proposed_radius);
        proposed_radius = (std::max)(b.zmax - b.zmin,
                                     proposed_radius);
        proposed_radius /= 100;
    }
    double r = QInputDialog::getDouble(NULL,
                                       tr("Display corners with new radius..."),
                                       tr("Radius:"),
                                       proposed_radius, // value
                                       0.,          // min
                                       2147483647., // max
                                       10,          // decimals
                                       &ok);
    if(ok) {
        change_corner_radii(r);
    }
}

void Scene_polylines_item::change_corner_radii(double r) {
    if(r >= 0) {
        d->spheres_drawn_radius = r;
        d->draw_extremities = (r > 0);
        this->changed();
        emit itemChanged();
    }
}

void Scene_polylines_item::split_at_sharp_angles()
{
    typedef Polylines_container Bare_polyline_container;
    typedef Polyline Bare_polyline;
    Polylines_container& bare_polylines = polylines;

    int counter = 0;
    for(Bare_polyline_container::iterator
        bare_polyline_it = bare_polylines.begin();
        bare_polyline_it != bare_polylines.end(); // the end changes
        // during the loop
        ++counter /* bare_polyline_it is incremented in the loop */)
    {
        Bare_polyline_container::iterator current_polyline_it =
                bare_polyline_it;
        Bare_polyline& bare_polyline = *bare_polyline_it;
        Bare_polyline::iterator it = boost::next(bare_polyline.begin());

        if(boost::next(bare_polyline.begin()) == bare_polyline.end())
        {
            std::cerr << "WARNING: Isolated point in polylines\n";
            bare_polyline_it = bare_polylines.erase(bare_polyline_it);
            continue;
        }
        else
            ++bare_polyline_it;
        if(it != bare_polyline.end()) {
            for(; it != boost::prior(bare_polyline.end()); ++it) {
                const Point_3 pv = *it;
                const Point_3 pa = *boost::prior(it);
                const Point_3 pb = *boost::next(it);
                const K::Vector_3 av = pv - pa;
                const K::Vector_3 bv = pv - pb;
                const K::FT sc_prod = av * bv;
                if( sc_prod >= 0 ||
                        (sc_prod < 0 &&
                         CGAL::square(sc_prod) < (av * av) * (bv * bv) / 4 ) )
                {
#ifdef PROTECTION_DEBUG
                    std::cerr << "Split polyline (small angle) "
                              <<  std::acos(sqrt(CGAL::square(sc_prod) /
                                                 ((av*av) * (bv*bv)))) * 180 /CGAL_PI
                               << " degres\n";
#endif
                    Bare_polyline new_polyline;
                    std::copy(it, bare_polyline.end(),
                              std::back_inserter(new_polyline));

                    if(*bare_polyline.begin() == *bare_polyline.rbegin()) {
                        // if the polyline is a cycle, test if its beginning is a sharp
                        // angle...
                        const Point_3 pv = *bare_polyline.begin();
                        const Point_3 pa = *boost::prior(boost::prior(bare_polyline.end()));
                        const Point_3 pb = *boost::next(bare_polyline.begin());
                        const K::Vector_3 av = pv - pa;
                        const K::Vector_3 bv = pv - pb;
                        const K::FT sc_prod = av * bv;
                        if( sc_prod >= 0 ||
                                (sc_prod < 0 &&
                                 CGAL::square(sc_prod) < (av * av) * (bv * bv) / 4 ) )
                        {
                            // if its beginning is a sharp angle, then split
                            bare_polyline.erase(boost::next(it), bare_polyline.end());
                        }
                        else {
                            // ...if not, modifies its beginning
                            std::copy(boost::next(bare_polyline.begin()),
                                      boost::next(it),
                                      std::back_inserter(new_polyline));
                            bare_polylines.erase(current_polyline_it);
                        }
                    }
                    else {
                        bare_polyline.erase(boost::next(it), bare_polyline.end());
                    }
                    bare_polylines.push_back(new_polyline);
                    break;
                }
            }
        }
    }
    emit itemChanged();
}

void
Scene_polylines_item::merge(Scene_polylines_item* other_item) {
    if(other_item == 0) return;
    std::copy(other_item->polylines.begin(),
              other_item->polylines.end(),
              std::back_inserter(polylines));
    QVariant other_metadata_variant = other_item->property("polylines metadata");
    if(other_metadata_variant.type() == QVariant::StringList)
    {
        QStringList metadata = property("polylines metadata").toStringList();
        metadata.append(other_metadata_variant.toStringList());
        setProperty("polylines metadata", metadata);
    }
    changed();
}

#include "Scene_polylines_item.moc"
