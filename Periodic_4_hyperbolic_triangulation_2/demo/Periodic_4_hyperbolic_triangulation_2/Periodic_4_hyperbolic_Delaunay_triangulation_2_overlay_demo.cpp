#include <fstream>

// CGAL headers
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_2.h> 
#include <CGAL/Periodic_4_hyperbolic_Delaunay_triangulation_traits_2.h> 
#include <CGAL/point_generators_2.h>
#include <CGAL/Hyperbolic_octagon_translation_matrix.h>

// unique words
#include <CGAL/Qt/HyperbolicPainterOstream.h>
#include <CGAL/Qt/utility.h>

// Qt headers
#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QGraphicsEllipseItem>

// for filtering
#include <set>
#include <string>

// GraphicsView items and event filters (input classes)
#include <CGAL/Qt/TriangulationCircumcircle.h>
#include <CGAL/Qt/TriangulationPointInputAndConflictZone.h>
#include <CGAL/Qt/TriangulationGraphicsItemWithColorInfoOverlay.h>     // Visualise color
#include <CGAL/Qt/DemosMainWindow.h>


#include <CGAL/CORE_Expr.h>
#include <CGAL/Cartesian.h>

// the two base classes
#include "ui_Periodic_4_hyperbolic_Delaunay_triangulation_2.h"



typedef CORE::Expr                                                              NT;
typedef CGAL::Cartesian<NT>                                                     Kernel;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_traits_2<Kernel>     Traits;
typedef CGAL::Periodic_4_hyperbolic_Delaunay_triangulation_2<Traits>            Triangulation;
typedef Hyperbolic_octagon_translation_matrix<Traits>                           Octagon_matrix;
typedef Kernel::Point_2                                                         Point;
typedef Triangulation::Vertex_handle                                            Vertex_handle;
typedef Traits::Side_of_fundamental_octagon                                     Side_of_fundamental_octagon;



class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Periodic_4_hyperbolic_Delaunay_triangulation_2
{
  Q_OBJECT
  
private:  


  int                                                           cidx;
  std::vector<int>                                              ccol;
  bool                                                          dummy_mode;

  Triangulation                                                 dt;
  QGraphicsEllipseItem                                        * disk;
  QGraphicsScene                                                scene;  

  CGAL::Qt::TriangulationGraphicsItem<Triangulation>          * dgi;
  
  CGAL::Qt::TriangulationPointInputAndConflictZone<Triangulation>  * pi;
  CGAL::Qt::TriangulationCircumcircle<Triangulation>               * tcc;
public:
  MainWindow();

public slots:

  void processInput(CGAL::Object o);
  
  void on_actionCircumcenter_toggled(bool checked);

  void on_actionShowTriangulation_toggled(bool checked);

  void on_actionInsertPoint_toggled(bool checked);
  
  void on_actionInsertRandomPoints_triggered();

  void on_actionLoadPoints_triggered();

  void on_actionClear_triggered();

  void on_actionRecenter_triggered();


signals:
  void changed();
};


MainWindow::MainWindow()
  : DemosMainWindow(), dt(Traits())
{

  dt.insert_dummy_points();

  cidx = 0;
  for (int i = 0; i < 14; i++)
    ccol.push_back(i);
  
  setupUi(this);

  this->graphicsView->setAcceptDrops(false);
  
  // Add Poincaré disk
  qreal origin_x = 0, origin_y = 0, radius = 1, diameter = 2*radius;
  qreal left_top_corner_x = origin_x - radius;
  qreal left_top_corner_y = origin_y - radius;
  qreal width = diameter, height = diameter;
  
  // set background
  qreal eps = 0.01;
  QGraphicsRectItem* rect = new QGraphicsRectItem(left_top_corner_x - eps, left_top_corner_y - eps, width + 2*eps, height + 2*eps);
  rect->setPen(Qt::NoPen);
  rect->setBrush(Qt::white);
  scene.addItem(rect);
  
  // set disk
  disk = new QGraphicsEllipseItem(left_top_corner_x, left_top_corner_y, width, height);
  QPen pen;  // creates a default pen
  pen.setWidth(0);
  //pen.setBrush(Qt::black);
  pen.setBrush(Qt::black);
  disk->setPen(pen);

  scene.addItem(disk);
  
  // Add a GraphicItem for the Triangulation triangulation
  dgi = new CGAL::Qt::TriangulationGraphicsItem<Triangulation>(&dt);

  QObject::connect(this, SIGNAL(changed()),
		   dgi, SLOT(modelChanged()));

  dgi->setVerticesPen(QPen(Qt::red, 3, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  dgi->setEdgesPen(QPen(QColor(200, 200, 0), 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  scene.addItem(dgi);

  // Setup input handlers. They get events before the scene gets them
  // and the input they generate is passed to the triangulation with 
  // the signal/slot mechanism    
  pi = new CGAL::Qt::TriangulationPointInputAndConflictZone<Triangulation>(&scene, &dt, this );

  QObject::connect(pi, SIGNAL(generate(CGAL::Object)),
		   this, SLOT(processInput(CGAL::Object)));

  tcc = new CGAL::Qt::TriangulationCircumcircle<Triangulation>(&scene, &dt, this);
  tcc->setPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));

  // 
  // Manual handling of actions
  //

  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   this, SLOT(close()));

  // We put mutually exclusive actions in an QActionGroup
  QActionGroup* ag = new QActionGroup(this);
  ag->addAction(this->actionInsertPoint);
  ag->addAction(this->actionCircumcenter);

  // Check two actions 
  this->actionInsertPoint->setChecked(true);
  this->actionShowTriangulation->setChecked(true);

  // //
  // // Setup the scene and the view
  // //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(left_top_corner_x, left_top_corner_y, width, height);
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);
  this->graphicsView->shear(230, 230);
  this->graphicsView->rotate(90);

  // // The navigation adds zooming and translation functionality to the
  // // QGraphicsView
  this->addNavigation(this->graphicsView);
  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutDemo(":/cgal/help/about_Triangulation_triangulation_2.html");
  this->addAboutCGAL();

}


void
MainWindow::processInput(CGAL::Object o)
{

  Point p;
  if(CGAL::assign(p, o)){
    Vertex_handle v = dt.insert(p);
  }
  emit(changed());

}




/* 
 *  Qt Automatic Connections
 *  http://doc.trolltech.com/4.4/designer-using-a-component.html#automatic-connections
 * 
 *  setupUi(this) generates connections to the slots named
 *  "on_<action_name>_<signal_name>"
 */
void
MainWindow::on_actionInsertPoint_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(pi);
  } else {
    scene.removeEventFilter(pi);
  }
}



void
MainWindow::on_actionCircumcenter_toggled(bool checked)
{
  if(checked){
    scene.installEventFilter(tcc);
    tcc->show();
  } else {  
    scene.removeEventFilter(tcc);
    tcc->hide();
  }
}



void
MainWindow::on_actionShowTriangulation_toggled(bool checked)
{
  dgi->setVisibleEdges(checked);
}



void
MainWindow::on_actionClear_triggered()
{
  dt.clear();
  dt.insert_dummy_points();
  emit(changed());
}


void
MainWindow::on_actionInsertRandomPoints_triggered()
{
  bool ok = false;
  const int number_of_points = 
    QInputDialog::getInt(this, 
                        tr("Number of random points"),
                        tr("Enter number of random points"),
           100,
           0,
           std::numeric_limits<int>::max(),
           1,
           &ok);

  if(!ok) {
    return;
  }

  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);


  typedef CGAL::Creator_uniform_2<double, Point> Creator;
  CGAL::Random_points_in_disc_2<Point, Creator> g( 1.0 );
  Side_of_fundamental_octagon pred;

  int cnt = 0;
  do {
    Point pt = *g;
    ++g;
    if (pred(pt) != CGAL::ON_UNBOUNDED_SIDE) {
      processInput(make_object(pt));
      cnt++;
    }
  } while (cnt < number_of_points);
  
  QApplication::restoreOverrideCursor();
  emit(changed());
}


void
MainWindow::on_actionLoadPoints_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open Points file"),
						  ".");
  if(! fileName.isEmpty()){
    //open(fileName);
  }
}



void
MainWindow::on_actionRecenter_triggered()
{
  this->graphicsView->setSceneRect(dgi->boundingRect());
  this->graphicsView->fitInView(dgi->boundingRect(), Qt::KeepAspectRatio);  
}


#include "Periodic_4_hyperbolic_Delaunay_triangulation_2_overlay_demo.moc"

int main(int argc, char **argv)
{

  QApplication app(argc, argv);
  app.setOrganizationDomain("geometryfactory.com");
  app.setOrganizationName("GeometryFactory");
  app.setApplicationName("Periodic_4_hyperbolic_Delaunay_triangulation_2_overlay demo");
  MainWindow mainWindow;
  mainWindow.show();
  QStringList args = app.arguments();
  args.removeAt(0);

  return app.exec();
}
