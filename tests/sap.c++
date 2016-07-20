#include "algorithm/SAP.h"
#include "debug/debug.h"


int main()
{
    DebugBegin();
    SAP<int> sap;
    AABB<2> bound1(1, 5, 1, 5);
    AABB<2> bound2(2, 10, 6, 11);
    AABB<2> bound3(6, 9, -2, 3);
    AABB<2> bound4(-10, -5, 0, 4);
    int lala;
    sap.add_box(bound1, &lala);
    sap.add_box(bound2, &lala);
    sap.add_box(bound3, &lala);
    sap.add_box(bound4, &lala);
    sap.update_box(1, AABB<2>(2, 10, 4, 11));
    sap.print();
    sap.assert_validity();
    sap.pairs.print();
    
}
