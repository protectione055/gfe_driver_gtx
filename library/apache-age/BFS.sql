SELECT * FROM ag_catalog.create_graph('my_graph');

-- 创建节点
SELECT * FROM cypher('my_graph', $$
    CREATE (a:Node {extid: 1, name: 'Node 1'}),
           (b:Node {extid: 2, name: 'Node 2'}),
           (c:Node {extid: 3, name: 'Node 3'}),
           (d:Node {extid: 4, name: 'Node 4'}),
           (e:Node {extid: 5, name: 'Node 5'})
$$) AS (result agtype);

-- 创建边（无向图需要两条方向相反的边）
-- 创建双向边
SELECT * FROM cypher('my_graph', $$
    MATCH (a:Node {extid: 1}), (b:Node {extid: 2})
    CREATE (a)-[:CONNECTED_TO]->(b),
           (b)-[:CONNECTED_TO]->(a)
$$) AS (result agtype);

SELECT * FROM cypher('my_graph', $$
    MATCH (b:Node {extid: 1}), (c:Node {extid: 3})
    CREATE (b)-[:CONNECTED_TO]->(c),
           (c)-[:CONNECTED_TO]->(b)
$$) AS (result agtype);

SELECT * FROM cypher('my_graph', $$
    MATCH (c:Node {extid: 2}), (d:Node {extid: 4})
    CREATE (c)-[:CONNECTED_TO]->(d),
           (d)-[:CONNECTED_TO]->(c)
$$) AS (result agtype);

SELECT * FROM cypher('my_graph', $$
    MATCH (d:Node {extid: 4}), (e:Node {extid: 5})
    CREATE (d)-[:CONNECTED_TO]->(e),
           (e)-[:CONNECTED_TO]->(d)
$$) AS (result agtype);


-- 执行BFS
WITH RECURSIVE bfs AS (
    SELECT id, extid, 1 AS level, ARRAY[id] AS visited
    FROM cypher('my_graph', $$
        MATCH (n:Node {extid: 1})
        RETURN id(n) AS id, n.extid AS extid
    $$) AS (id agtype, extid agtype)
    
    UNION ALL
    
    SELECT e.end_id AS id, e.end_extid AS extid, p.level + 1 AS level, p.visited || e.end_id AS visited
    FROM bfs AS p
    JOIN cypher('my_graph', $$
        MATCH (n:Node)-[r]->(m:Node)
        RETURN id(n) AS start_id, id(m) AS end_id, m.extid AS end_extid
    $$) AS e(start_id agtype, end_id agtype, end_extid agtype) ON e.start_id = p.id
    WHERE e.end_id<>ALL(p.visited)
)
SELECT * FROM bfs;

SELECT * FROM cypher('my_graph', $$
    MATCH (n:Node {id: 1})-[:CONNECTED_TO*]->(m)
    RETURN m
$$) AS (m agtype);
