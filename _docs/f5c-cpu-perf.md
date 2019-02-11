---
title: F5C multi-core CPU performance
chart: true
---

## Runtime for processing (excluding I/O)

The dataset was copied to `/dev/shm1` (tmpfs file system that resides in RAM) so that the disk I/O time is excluded[^1]. Runtime measurements were taken on a server (72 Intel Xeon Gold 6154 threads and 376 GB RAM).

Performance of F5C over Nanopolish call-methylation v0.9.0 (Nanopolish version which the F5C development was based on) is as follows :

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
      "text": "Execution time with respect to number of threads (data on RAM)"
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

Note that F5C was run with the following additional options to tally with the default hard coded parameters in Nanopolish.
- min-mapq: 0
- batch size: 512

We observe an approximately 1.5x to 2x performance improvement over Nanopolish when using f5c. Based on what was learnt during F5C development, a number of improvements were applied to the original Nanopolish repository : [#486](https://github.com/jts/nanopolish/pull/486),[#350](https://github.com/jts/nanopolish/pull/350),[#402](https://github.com/jts/nanopolish/pull/402),[#327](https://github.com/jts/nanopolish/pull/327),[#285](https://github.com/jts/nanopolish/issues/285), [364](https://github.com/jts/nanopolish/issues/364). Following is a comparison with Nanopolish  v0.10.2 with those improvements added.


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
        "label": "nanopolish(v0.10.2) call-methylation",
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
      "text": "Execution time with respect to number of threads (data on RAM)"
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

## Total runtime (with disk I/O)

Instead of loading data from the RAM, now we directly loaded from the hard disk drives. The same above server used for taking measurements had a RAID 6 setup of 12 mechanical hard disk drives. We cleared the operating system file system cache before each measurement using `sync; echo 3 | tee /proc/sys/vm/drop_caches`.

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
        "label": "nanpolish(v0.10.2) call-methylation",
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

It is observed that the performance gain for a large number of threads is limited. It is mostly due to disk I/O becomes performance bottleneck, where you would not see any performance improvement because loading the files into memory takes more time than the actual processing. See [HDF5 Performance](/docs/hdf5-performance) about the impact of the I/O.

[^1]: even when loading from tmpfs (RAM), data loading time of F5C is higher (almost twice) the processing time as shown below.
```
[meth_main] Data loading time: 43.642 sec
[meth_main]     - bam load time: 1.662 sec
[meth_main]     - fasta load time: 11.611 sec
[meth_main]     - fast5 load time: 30.305 sec
[meth_main]         - fast5 open time: 1.775 sec
[meth_main]         - fast5 read time: 25.027 sec
[meth_main] Data processing time: 24.404 sec
```
This might be due to the overheads of the library calls for accessing files. As a result, the multi-threaded performance of F5C plateaus after around 32 threads.
