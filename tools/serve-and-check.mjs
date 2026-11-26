/**
 * Serve dist/ on localhost and run the browser checks against it.
 *
 * Bound to 127.0.0.1 deliberately: the built site is served only to this machine,
 * never exposed on a public interface.
 */
import { spawn } from "node:child_process";
import { createServer } from "node:http";
import { readFile, stat } from "node:fs/promises";
import { extname, join, normalize } from "node:path";

const PORT = Number(process.env.PORT ?? 8099);
const ROOT = new URL("../dist/", import.meta.url).pathname;
const TYPES = {
  ".html": "text/html", ".js": "text/javascript", ".css": "text/css",
  ".json": "application/json", ".bin": "application/octet-stream",
  ".idx": "application/octet-stream", ".dir": "application/octet-stream",
};

const server = createServer(async (req, res) => {
  try {
    const path = decodeURIComponent(new URL(req.url, "http://x").pathname);
    const file = join(ROOT, normalize(path).replace(/^(\.\.[/\\])+/, ""));
    const target = (await stat(file)).isDirectory() ? join(file, "index.html") : file;
    const body = await readFile(target);
    const range = req.headers.range;
    const match = range && /bytes=(\d+)-(\d+)?/.exec(range);
    const type = TYPES[extname(target)] ?? "application/octet-stream";

    if (match) {
      const start = Number(match[1]);
      const end = match[2] ? Number(match[2]) : body.length - 1;
      res.writeHead(206, {
        "Content-Type": type,
        "Content-Range": `bytes ${start}-${end}/${body.length}`,
        "Accept-Ranges": "bytes",
      });
      res.end(body.subarray(start, end + 1));
      return;
    }
    res.writeHead(200, { "Content-Type": type, "Accept-Ranges": "bytes" });
    res.end(body);
  } catch {
    res.writeHead(404).end("not found");
  }
});

server.listen(PORT, "127.0.0.1", () => {
  const shots = process.env.SHOT_DIR ?? "/tmp/splashrna-shots";
  const child = spawn(
    process.execPath,
    [new URL("browser-check.mjs", import.meta.url).pathname, `http://127.0.0.1:${PORT}/`, shots],
    { stdio: "inherit" },
  );
  child.on("exit", (code) => {
    server.close();
    process.exit(code ?? 1);
  });
});
