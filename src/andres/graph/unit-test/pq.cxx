#include <stdexcept>
#include <iostream>
#include "andres/graph/pq.hxx"
#include "andres/graph/runtime_check.hxx"





void testMinQueue(){
    const float tol=0.001f;
    {
        typedef graph::ChangeablePriorityQueue<float,std::less<float> > MinQueueType;
        auto q = MinQueueType(4);

        auto b = q.empty();
        // CURRENT VALUES
        //-----------------
        // 0 => NONE
        // 1 => NONE
        // 2 => NONE
        // 3 => NONE
        GRAPH_TEST(q.empty());
        GRAPH_TEST_OP(q.size(),==,0);
        GRAPH_TEST(!q.contains(0));
        GRAPH_TEST(!q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST(!q.contains(3));



        q.push(0,3.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => NONE
        // 2 => NONE
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,1);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST(!q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(1),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(2),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),3.0, tol);
        GRAPH_TEST_OP(q.top(),==,0);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),3.0,tol);


        q.push(2,2.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => NONE
        // 2 => 2.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,2);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST(!q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(1),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),3.0, tol);
        GRAPH_TEST_OP(q.top(),==,2);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),2.0,tol);


        q.push(1,3.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 3.0
        // 2 => 2.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,3);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),3.0, tol);
        GRAPH_TEST_OP(q.top(),==,2);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),2.0,tol);


        q.push(3,0.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 3.0
        // 2 => 2.0
        // 3 => 0.0
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,4);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST( q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(3),0.0, tol);
        GRAPH_TEST_OP(q.top(),==,3);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),0.0,tol);

        q.pop();
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 3.0
        // 2 => 2.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,3);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),0.0, 0.01);
        GRAPH_TEST_OP(q.top(),==,2);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),2.0,tol);


        q.push(1,1.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 1.0
        // 2 => 2.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,3);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),1.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),0.0, 0.01);
        GRAPH_TEST_OP(q.top(),==,1);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),1.0,tol);

        q.push(1,4.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 4.0
        // 2 => 2.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,3);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),0.0, 0.01);
        GRAPH_TEST_OP(q.top(),==,2);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),2.0,tol);


        q.push(0,1.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 1.0
        // 1 => 4.0
        // 2 => 2.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,3);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),1.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),0.0, 0.01);
        GRAPH_TEST_OP(q.top(),==,0);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),1.0,tol);


        q.pop();
        // CURRENT VALUES
        //-----------------
        // 0 => NONE
        // 1 => 4.0
        // 2 => 2.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,2);
        GRAPH_TEST(!q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        //GRAPH_CHECK_EQ_TOL(q.priority(0),1.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),0.0, 0.01);
        GRAPH_TEST_OP(q.top(),==,2);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),2.0,tol);


        q.pop();
        // CURRENT VALUES
        //-----------------
        // 0 => NONE
        // 1 => 4.0
        // 2 => NONE
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,1);
        GRAPH_TEST(!q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST(!q.contains(3));
        //GRAPH_CHECK_EQ_TOL(q.priority(0),1.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),0.0, 0.01);
        GRAPH_TEST_OP(q.top(),==,1);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),4.0,tol);


        q.pop();
        // CURRENT VALUES
        //-----------------
        // 0 => NONE
        // 1 => NONE
        // 2 => NONE
        // 3 => NONE
        GRAPH_TEST(q.empty());
        GRAPH_TEST_OP(q.size(),==,0);
        GRAPH_TEST(!q.contains(0));
        GRAPH_TEST(!q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST(!q.contains(3));
        //GRAPH_CHECK_EQ_TOL(q.priority(0),1.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),0.0, 0.01);
        //GRAPH_TEST_OP(q.top(),==,1);
        //GRAPH_CHECK_EQ_TOL(q.topPriority(),4.0,tol);

    }

}

void testMaxQueue(){
    const float tol=0.001f;
    {
        typedef graph::ChangeablePriorityQueue<float,std::greater<float> > MaxQueueType;
        MaxQueueType q(4);

        // CURRENT VALUES
        //-----------------
        // 0 => NONE
        // 1 => NONE
        // 2 => NONE
        // 3 => NONE
        GRAPH_TEST(q.empty());
        GRAPH_TEST_OP(q.size(),==,0);
        GRAPH_TEST(!q.contains(0));
        GRAPH_TEST(!q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST(!q.contains(3));



        q.push(0,3.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => NONE
        // 2 => NONE
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,1);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST(!q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(1),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(2),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),3.0, tol);
        GRAPH_TEST_OP(q.top(),==,0);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),3.0,tol);


        q.push(2,2.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => NONE
        // 2 => 2.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,2);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST(!q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(1),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),3.0, tol);
        GRAPH_TEST_OP(q.top(),==,0);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),3.0,tol);


        q.push(1,4.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 4.0
        // 2 => 2.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,3);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),3.0, tol);
        GRAPH_TEST_OP(q.top(),==,1);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),4.0,tol);


        q.push(3,5.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 4.0
        // 2 => 2.0
        // 3 => 5.0
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,4);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST( q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(3),5.0, tol);
        GRAPH_TEST_OP(q.top(),==,3);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),5.0,tol);


        q.push(3,2.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 4.0
        // 2 => 2.0
        // 3 => 2.0
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,4);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST( q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(3),2.0, tol);
        GRAPH_TEST_OP(q.top(),==,1);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),4.0,tol);


        q.push(1,0.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 0.0
        // 2 => 2.0
        // 3 => 2.0
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,4);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST( q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),0.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),2.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(3),2.0, tol);
        GRAPH_TEST_OP(q.top(),==,0);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),3.0,tol);

        q.push(2,5.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 0.0
        // 2 => 5.0
        // 3 => 2.0
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,4);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST( q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),0.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),5.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(3),2.0, tol);
        GRAPH_TEST_OP(q.top(),==,2);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),5.0,tol);

        q.push(3,1.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 0.0
        // 2 => 5.0
        // 3 => 1.0
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,4);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST( q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),0.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),5.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(3),1.0, tol);
        GRAPH_TEST_OP(q.top(),==,2);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),5.0,tol);


        q.pop();
        // CURRENT VALUES
        //-----------------
        // 0 => 3.0
        // 1 => 0.0
        // 2 => NONE
        // 3 => 1.0
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,3);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST( q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),0.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(2),5.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(3),1.0, tol);
        GRAPH_TEST_OP(q.top(),==,0);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),3.0,tol);



        q.pop();
        // CURRENT VALUES
        //-----------------
        // 0 => NONE
        // 1 => 0.0
        // 2 => NONE
        // 3 => 1.0
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,2);
        GRAPH_TEST(!q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST( q.contains(3));
        //GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),0.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(2),5.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(3),1.0, tol);
        GRAPH_TEST_OP(q.top(),==,3);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),1.0,tol);


        q.pop();
        // CURRENT VALUES
        //-----------------
        // 0 => NONE
        // 1 => 0.0
        // 2 => NONE
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,1);
        GRAPH_TEST(!q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST(!q.contains(3));
        //GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),0.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(2),5.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),1.0, tol);
        GRAPH_TEST_OP(q.top(),==,1);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),0.0,tol);

        q.pop();
        // CURRENT VALUES
        //-----------------
        // 0 => NONE
        // 1 => NONE
        // 2 => NONE
        // 3 => NONE
        GRAPH_TEST(q.empty());
        GRAPH_TEST_OP(q.size(),==,0);
        GRAPH_TEST(!q.contains(0));
        GRAPH_TEST(!q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST(!q.contains(3));
        //GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(1),0.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(2),5.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),1.0, tol);
        //GRAPH_TEST_OP(q.top(),==,1);
        //GRAPH_CHECK_EQ_TOL(q.topPriority(),0.0,tol);

        q.push(2,1.0);
        // CURRENT VALUES
        //-----------------
        // 0 => NONE
        // 1 => NONE
        // 2 => 1.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,1);
        GRAPH_TEST(!q.contains(0));
        GRAPH_TEST(!q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        //GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(1),0.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),1.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),1.0, tol);
        GRAPH_TEST_OP(q.top(),==,2);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),1.0,tol);


        q.push(2,3.0);
        // CURRENT VALUES
        //-----------------
        // 0 => NONE
        // 1 => NONE
        // 2 => 3.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,1);
        GRAPH_TEST(!q.contains(0));
        GRAPH_TEST(!q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        //GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(1),0.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),1.0, tol);
        GRAPH_TEST_OP(q.top(),==,2);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),3.0,tol);


        q.push(1,4.0);
        // CURRENT VALUES
        //-----------------
        // 0 => NONE
        // 1 => 4.0
        // 2 => 3.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,2);
        GRAPH_TEST(!q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        //GRAPH_CHECK_EQ_TOL(q.priority(0),3.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),1.0, tol);
        GRAPH_TEST_OP(q.top(),==,1);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),4.0,tol);


        q.push(0,2.0);
        // CURRENT VALUES
        //-----------------
        // 0 => 2.0
        // 1 => 4.0
        // 2 => 3.0
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,3);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST( q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),2.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(2),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),1.0, tol);
        GRAPH_TEST_OP(q.top(),==,1);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),4.0,tol);


        q.deleteItem(2);
        // CURRENT VALUES
        //-----------------
        // 0 => 2.0
        // 1 => 4.0
        // 2 => NONE
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,2);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST( q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),2.0, tol);
        GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(2),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),1.0, tol);
        GRAPH_TEST_OP(q.top(),==,1);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),4.0,tol);


        q.deleteItem(1);
        // CURRENT VALUES
        //-----------------
        // 0 => 2.0
        // 1 => 4.0
        // 2 => NONE
        // 3 => NONE
        GRAPH_TEST(!q.empty());
        GRAPH_TEST_OP(q.size(),==,1);
        GRAPH_TEST( q.contains(0));
        GRAPH_TEST(!q.contains(1));
        GRAPH_TEST(!q.contains(2));
        GRAPH_TEST(!q.contains(3));
        GRAPH_CHECK_EQ_TOL(q.priority(0),2.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(1),4.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(2),3.0, tol);
        //GRAPH_CHECK_EQ_TOL(q.priority(3),1.0, tol);
        GRAPH_TEST_OP(q.top(),==,0);
        GRAPH_CHECK_EQ_TOL(q.topPriority(),2.0,tol);

    }

}

int main() {
    testMinQueue();
    testMaxQueue();
    return 0;
}
