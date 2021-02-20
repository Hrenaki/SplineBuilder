using Microsoft.VisualStudio.TestTools.UnitTesting;
using System;
using SplineBuilder;
using System.Collections.Generic;

namespace SplineTest
{
    [TestClass]
    public class UnitTest1
    {
        [TestMethod]
        public void TestMethod1()
        {
            List<float[]> splinePoints = new List<float[]>() { new float[] { -1.0f, 0.0f }, new float[]{ 0.0f, 0.0f},
                new float[]{ 1.0f, 0.0f} };
            Spline s = new Spline(splinePoints);
            s.CreateSpline();
        }
    }
}
