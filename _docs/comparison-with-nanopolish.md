---
title: Comparison with Nanopolish
chart: true
---

The disk cache could have significant impact on the execution time of both f5c
and nanopolish. Running without disk cache (for example, first f5c run) could be
a performance bottleneck when number of threads is large enough, where you would
not see any performance improvement because loading the files into memory takes
more time than the actual processing.

```chart
{
  "type": "line",
  "data": {
    "labels": [
      "8",
      "16",
      "24",
      "32",
      "40",
      "48",
      "56",
      "64",
      "72"
    ],
    "datasets": [
      {
        "label": "f5c",
        "fill": false,
        "lineTension": 0.1,
        "backgroundColor": "rgba(75,192,192,0.4)",
        "borderColor": "rgba(75,192,192,1)",
        "borderCapStyle": "butt",
        "borderDash": [],
        "borderDashOffset": 0,
        "borderJoinStyle": "miter",
        "pointBorderColor": "rgba(75,192,192,1)",
        "pointBackgroundColor": "#fff",
        "pointBorderWidth": 1,
        "pointHoverRadius": 5,
        "pointHoverBackgroundColor": "rgba(75,192,192,1)",
        "pointHoverBorderColor": "rgba(220,220,220,1)",
        "pointHoverBorderWidth": 2,
        "pointRadius": 1,
        "pointHitRadius": 10,
        "data": [
          163.09,
          97.02,
          74.64,
          70.24,
          65.78,
          67.08,
          64.85,
          65.5,
          66.19
        ],
        "spanGaps": false
      },
      {
        "label": "nanpolish call-methylation",
        "fill": false,
        "lineTension": 0.1,
        "backgroundColor": "rgba(255, 99, 132, 0.2)",
        "borderColor": "rgba(255,99,132,1)",
        "borderCapStyle": "butt",
        "borderDash": [],
        "borderDashOffset": 0,
        "borderJoinStyle": "miter",
        "pointBorderColor": "#fff",
        "pointBackgroundColor": "rgba(255,99,132,1)",
        "pointBorderWidth": 1,
        "pointHoverRadius": 5,
        "pointHoverBackgroundColor": "#fff",
        "pointHoverBorderColor": "rgba(255,99,132,1)",
        "pointHoverBorderWidth": 2,
        "pointRadius": 1,
        "pointHitRadius": 10,
        "data": [
          198,
          111.78,
          86.61,
          79.14,
          81,
          95.61,
          99.59,
          101.29,
          101.53
        ],
        "spanGaps": false
      }
    ]
  },
  "options": {
    "title": {
      "display": true,
      "text": "Execution time with respect to number of threads"
    },
    "scales": {
      "yAxes": [{
        "scaleLabel": {
          "display": true,
          "labelString": "time/s"
        }
      }],
      "xAxes": [{
        "scaleLabel": {
          "display": true,
          "labelString": "threads"
        }
      }]
    }
  }
}
```

```chart
{
  "type": "line",
  "data": {
    "labels": [
      "8",
      "16",
      "24",
      "32",
      "40",
      "48",
      "56",
      "64",
      "72"
    ],
    "datasets": [
      {
        "label": "f5c",
        "fill": false,
        "lineTension": 0.1,
        "backgroundColor": "rgba(75,192,192,0.4)",
        "borderColor": "rgba(75,192,192,1)",
        "borderCapStyle": "butt",
        "borderDash": [],
        "borderDashOffset": 0,
        "borderJoinStyle": "miter",
        "pointBorderColor": "rgba(75,192,192,1)",
        "pointBackgroundColor": "#fff",
        "pointBorderWidth": 1,
        "pointHoverRadius": 5,
        "pointHoverBackgroundColor": "rgba(75,192,192,1)",
        "pointHoverBorderColor": "rgba(220,220,220,1)",
        "pointHoverBorderWidth": 2,
        "pointRadius": 1,
        "pointHitRadius": 10,
        "data": [
          172.23,
          161.02,
          161.93,
          163.3,
          164.45,
          164.49,
          161.69,
          164.33,
          161.32
        ],
        "spanGaps": false
      },
      {
        "label": "nanpolish call-methylation",
        "fill": false,
        "lineTension": 0.1,
        "backgroundColor": "rgba(255, 99, 132, 0.2)",
        "borderColor": "rgba(255,99,132,1)",
        "borderCapStyle": "butt",
        "borderDash": [],
        "borderDashOffset": 0,
        "borderJoinStyle": "miter",
        "pointBorderColor": "#fff",
        "pointBackgroundColor": "rgba(255,99,132,1)",
        "pointBorderWidth": 1,
        "pointHoverRadius": 5,
        "pointHoverBackgroundColor": "#fff",
        "pointHoverBorderColor": "rgba(255,99,132,1)",
        "pointHoverBorderWidth": 2,
        "pointRadius": 1,
        "pointHitRadius": 10,
        "data": [
          240.73,
          177.72,
          185.99,
          183.68,
          186.64,
          193.32,
          195.66,
          198.39,
          200.18
        ],
        "spanGaps": false
      }
    ]
  },
  "options": {
    "title": {
      "display": true,
      "text": "Execution time with respect to number of threads (without disk cache)"
    },
    "scales": {
      "yAxes": [{
        "scaleLabel": {
          "display": true,
          "labelString": "time/s"
        }
      }],
      "xAxes": [{
        "scaleLabel": {
          "display": true,
          "labelString": "threads"
        }
      }]
    }
  }
}
```

---

Performance can best be benchmarked by putting the data set in ramdisk, where
disk cache has no effect on the execution time of both f5c and nanopolish.

```chart
{
  "type": "line",
  "data": {
    "labels": [
      "8",
      "16",
      "24",
      "32",
      "40",
      "48",
      "56",
      "64",
      "72"
    ],
    "datasets": [
      {
        "label": "f5c",
        "fill": false,
        "lineTension": 0.1,
        "backgroundColor": "rgba(75,192,192,0.4)",
        "borderColor": "rgba(75,192,192,1)",
        "borderCapStyle": "butt",
        "borderDash": [],
        "borderDashOffset": 0,
        "borderJoinStyle": "miter",
        "pointBorderColor": "rgba(75,192,192,1)",
        "pointBackgroundColor": "#fff",
        "pointBorderWidth": 1,
        "pointHoverRadius": 5,
        "pointHoverBackgroundColor": "rgba(75,192,192,1)",
        "pointHoverBorderColor": "rgba(220,220,220,1)",
        "pointHoverBorderWidth": 2,
        "pointRadius": 1,
        "pointHitRadius": 10,
        "data": [
          199.55333333333337,
          162.99,
          162.67999999999998,
          164.84666666666666,
          163.24333333333334,
          164.60333333333332,
          162.71666666666667,
          163.0966666666667,
          168.4
        ],
        "spanGaps": false
      },
      {
        "label": "nanpolish(v0.9.0) call-methylation",
        "fill": false,
        "lineTension": 0.1,
        "backgroundColor": "rgba(255, 99, 132, 0.2)",
        "borderColor": "rgba(255,99,132,1)",
        "borderCapStyle": "butt",
        "borderDash": [],
        "borderDashOffset": 0,
        "borderJoinStyle": "miter",
        "pointBorderColor": "#fff",
        "pointBackgroundColor": "rgba(255,99,132,1)",
        "pointBorderWidth": 1,
        "pointHoverRadius": 5,
        "pointHoverBackgroundColor": "#fff",
        "pointHoverBorderColor": "rgba(255,99,132,1)",
        "pointHoverBorderWidth": 2,
        "pointRadius": 1,
        "pointHitRadius": 10,
        "data": [
          360.08,
          265.61,
          257.855,
          249.075,
          252.385,
          259.265,
          264.22,
          268.58,
          270.32
        ],
        "spanGaps": false
      }
    ]
  },
  "options": {
    "title": {
      "display": true,
      "text": "Execution time with respect to number of threads (data on HDD)"
    },
    "scales": {
      "yAxes": [{
        "scaleLabel": {
          "display": true,
          "labelString": "time/s"
        }
      }],
      "xAxes": [{
        "scaleLabel": {
          "display": true,
          "labelString": "threads"
        }
      }]
    }
  }
}
```

```chart
{
  "type": "line",
  "data": {
    "labels": [
      "8",
      "16",
      "24",
      "32",
      "40",
      "48",
      "56",
      "64",
      "72"
    ],
    "datasets": [
      {
        "label": "f5c",
        "fill": false,
        "lineTension": 0.1,
        "backgroundColor": "rgba(75,192,192,0.4)",
        "borderColor": "rgba(75,192,192,1)",
        "borderCapStyle": "butt",
        "borderDash": [],
        "borderDashOffset": 0,
        "borderJoinStyle": "miter",
        "pointBorderColor": "rgba(75,192,192,1)",
        "pointBackgroundColor": "#fff",
        "pointBorderWidth": 1,
        "pointHoverRadius": 5,
        "pointHoverBackgroundColor": "rgba(75,192,192,1)",
        "pointHoverBorderColor": "rgba(220,220,220,1)",
        "pointHoverBorderWidth": 2,
        "pointRadius": 1,
        "pointHitRadius": 10,
        "data": [
          162.54666666666668,
          97.38,
          75.58666666666666,
          73.31,
          69.76666666666667,
          70.58666666666666,
          68.69,
          69.41333333333334,
          69.64999999999999
        ],
        "spanGaps": false
      },
      {
        "label": "nanpolish(v0.9.0) call-methylation",
        "fill": false,
        "lineTension": 0.1,
        "backgroundColor": "rgba(255, 99, 132, 0.2)",
        "borderColor": "rgba(255,99,132,1)",
        "borderCapStyle": "butt",
        "borderDash": [],
        "borderDashOffset": 0,
        "borderJoinStyle": "miter",
        "pointBorderColor": "#fff",
        "pointBackgroundColor": "rgba(255,99,132,1)",
        "pointBorderWidth": 1,
        "pointHoverRadius": 5,
        "pointHoverBackgroundColor": "#fff",
        "pointHoverBorderColor": "rgba(255,99,132,1)",
        "pointHoverBorderWidth": 2,
        "pointRadius": 1,
        "pointHitRadius": 10,
        "data": [
          323.8833333333334,
          221.57333333333335,
          194.34,
          181.81333333333336,
          183.52333333333334,
          183.01,
          187.85999999999999,
          187.9333333333333,
          189.0233333333333
        ],
        "spanGaps": false
      }
    ]
  },
  "options": {
    "title": {
      "display": true,
      "text": "Execution time with respect to number of threads (data on ramdisk)"
    },
    "scales": {
      "yAxes": [{
        "scaleLabel": {
          "display": true,
          "labelString": "time/s"
        }
      }],
      "xAxes": [{
        "scaleLabel": {
          "display": true,
          "labelString": "threads"
        }
      }]
    }
  }
}
```

Now we can see that there is an approximately 1.5x to 2x performance improvement
when using f5c.
